/* ===================================================================
 * cu_fastq_stats.cc
 *  computes basic sequence stats for a FASTQ dataset using the GPU
 *  These include average A,T,C,G and N content and average quality
 *  per position in the sequence
 * Written by Assen Roguev, 2018
 * =================================================================== */
#include <stdio.h>
#include <cstdlib>
#include <cstdint>
#include <memory.h>
#include <errno.h>
#include <err.h>

// compile with -DCUDA_ERROR_CHECK toturn on error checking
#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line ) {
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err ) {
        fprintf( stderr, "%s:%i : %s\n", file, line, cudaGetErrorString( err ) );
        
        exit( -1 );
    }
#endif
    return;
}

inline void __cudaCheckError( const char *file, const int line ) {
#ifdef CUDA_ERROR_CHECK
    cudaError err = cudaGetLastError();
    if ( cudaSuccess != err ) {
        fprintf( stderr, "%s:%i : %s\n", cudaGetErrorString( err ) );
        exit( -1 );
    }

    // More careful checking. However, this will affect performance.
    // Comment away if needed.
    err = cudaDeviceSynchronize();
    if( cudaSuccess != err ) {
        fprintf( stderr, "%s:%i : %s\n", file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif
    return;
}

/* pos_stats is a cuda kernel to compuye sequence stats
 * Arguments:
 *  device char* to concatenated sequence data
 *  device float* to hold the stats
 *  size_t length of sequence read
 *  size_t chunk size 
 */
__global__
void pos_stats(char* str, float* pos_stats, size_t LEN, size_t CHUNK_SZ) {
    // define shared memory array
    extern __shared__ int spos_stats[];
    // set to 0
    if (threadIdx.x == 0) {
        for (size_t i = 0; i < 6*LEN; i++)
            spos_stats[i] = 0;
    }
    __syncthreads();
    
    // get the sequence index
    unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
    if (ix < CHUNK_SZ) {
        for (size_t i = 0; i < LEN; i++) {
            switch (str[2*LEN*ix+i]) {
                case 'A' :  case 'a':
                    atomicAdd(&spos_stats[6*i+0], 1);
                    break;
                case 'T' :  case 't':
                    atomicAdd(&spos_stats[6*i+1], 1);
                    break;
                case 'G' :  case 'g':
                    atomicAdd(&spos_stats[6*i+2], 1);
                    break;
                case 'C' :  case 'c':
                    atomicAdd(&spos_stats[6*i+3], 1);
                    break;
                case 'N' : case 'n':
                    atomicAdd(&spos_stats[6*i+4], 1);
                    break;
            }
        // per-position quality
        atomicAdd(&spos_stats[6*i+5], str[2*LEN*ix+i+LEN]);
        }
    }
    __syncthreads();
    
    // add to the global pos_stats structure
    if (threadIdx.x == 0) {
        for (size_t i = 0; i < 6*LEN; i++)
            atomicAdd(&pos_stats[i], spos_stats[i]);
    }
}

/* print_pos_stats prints the contents of the pos_stats array
 * Arguments:
 *  device float* containing the stats to be printed
 *  size_t length of the sequence read
 *  uint64_t total number of sequences
 */
void print_pos_stats(float* d_pos_stats, size_t LEN, uint64_t totalN) {
    const char charset[] = "ATGCNQ";
    
    // copy data to hist memory
    float* tmp_pos_stats = (float*)calloc(6*LEN, sizeof(float));
    CudaSafeCall( cudaMemcpy(tmp_pos_stats, d_pos_stats, 6*LEN*sizeof(float), cudaMemcpyDeviceToHost) );
    
    for (int b = 0; b < 6; b++) {
        printf("@%c", charset[b]);
        for (size_t i = 0; i < LEN; i++) {
            if (b < 5)
                // average nucleotide content
                printf("\t%0.4f",(float)tmp_pos_stats[6*i+b]/totalN);
            else {
                // check if q-stats are 0 at this position
                // can happen if the LEN is bigger than the actual read length
                if (tmp_pos_stats[6*i+b] == 0)
                     printf("\t%0.4f", (float)0);
                else
                    // average quality
                    printf("\t%0.4f",(float)(tmp_pos_stats[6*i+b] - 33*totalN)/totalN);
            }
        }
        printf("\n");
    }
    // cleanup
    free(tmp_pos_stats);
}

/* read_fastq attempts to read a fastq record from a file
 * Arguments:
 *  char* to hold the read data
 *  FILE* to an open file
 *  size_t for sequence read length
 *  size_t offset within the target string
 */
bool read_fastq(char* str, FILE* fp, const size_t L, size_t offset) {
    char* line = NULL;
    size_t len = 0;
    ssize_t read, rs;
    
    while (1) {
        // find entry point
        read = getline(&line, &len, fp);
        // EOF or other error
        if (read == -1)
            break;
        
        // entry point, line starts with '@'
        if (line[0] != '@') { 
            free(line); 
            continue;   // drop this line keep going
        }
            
        // read sequence string
        free(line);
        line = NULL;
        len = 0;
        read = getline(&line, &len, fp);
        // check for error and empty line
        if ((read != -1) && (read != 1)) {
            rs = read - 1;      // remove newline
            if (rs < L) 
                memcpy(str+2*L*offset,line,rs*sizeof(char)); 
            else
                memcpy(str+2*L*offset,line,L*sizeof(char));
        } else
            break;
            
        // read '+' string, line starts with '+'
        free(line);
        line = NULL;
        len = 0;
        read = getline(&line, &len, fp);
        // check for error and not starting with '+'
        if ((read != -1) && (line[0] != '+')) 
            break; 
            
        // read q-string
        free(line);
        line = NULL;
        len = 0;
        read = getline(&line, &len, fp);
        // check for error and different length from sequence
        if ((read != -1) && (read - 1 == rs)) {
            // newline eliminated already
            if (rs < L)
                memcpy(str+2*L*offset+L,line,rs*sizeof(char));
            else 
                memcpy(str+2*L*offset+L,line,L*sizeof(char));
            
            free(line);
            return true;    // success
        
        } else
            break;
    }
    
    // something went wrong
    free(line);
    return false;
}

/* main 
 * Arguments:
 *  filename
 *  sequence read length
 *  (optional) max number of sequences to process
 */
int main(int argc, char** argv) {
    // change this if needed (and change the kernel launch parameters below)
    size_t CHUNK_SZ = 1000000;
    
    if (argc < 3) {
        fprintf(stderr, "%s: Insufficient arguments\n", argv[0]);
        exit(1);
    }
    
    char* fname = argv[1];              // input filename
    size_t LEN = (size_t)atoi(argv[2]); // desired sequence length
    uint64_t maxSeq = 0;                // maximum number of sequences to process
    uint64_t totalN = 0;                // total number of sequences
    
    if (argc == 4) {
        maxSeq = (uint64_t)atoi(argv[3]);
    }
    
    // open file for reading
    FILE* fp = fopen(fname, "r");
    if (fp == NULL) { 
        err(1, "fopen: %s", fname);
        exit(2);
    }
    
    // allocate host and device memory
    char* h_str = (char*)calloc(2*LEN*CHUNK_SZ,sizeof(char));
    
    char* d_str;       // sequence
    CudaSafeCall( cudaMalloc((void**)&d_str, 2*CHUNK_SZ*LEN*sizeof(char)) );
    CudaSafeCall( cudaMemset(d_str, 0, 2*CHUNK_SZ*LEN*sizeof(char)) );


    float* d_pos_stats;  // pos stats    
    CudaSafeCall( cudaMalloc((void**)&d_pos_stats, 6*LEN*sizeof(float)) );
    CudaSafeCall( cudaMemset(d_pos_stats, 0, 6*LEN*sizeof(float)) );
    
    // main loop to go arbitrary number of sequences
    bool done = false;
    while(!done) {
        for (size_t i = 0; i < (size_t)CHUNK_SZ; i++) {
            if (read_fastq(h_str, fp, LEN, i)) {
                totalN++;
                if (totalN == maxSeq) {
                    done = true;
                    break;
                }
            
            } else { 
                done = true; 
                break;
            }
        }

        CudaSafeCall( cudaMemcpy(d_str, h_str, 2*CHUNK_SZ*LEN*sizeof(char), cudaMemcpyHostToDevice) );
    
        pos_stats<<<1024,1024,6*LEN*sizeof(float)>>>(d_str, d_pos_stats, LEN, CHUNK_SZ);
        CudaCheckError();
            
        // reset arrays
        memset(h_str, 0, 2*CHUNK_SZ*LEN*sizeof(char));
        CudaSafeCall( cudaMemset(d_str, 0, 2*CHUNK_SZ*LEN*sizeof(char)) );
    }
    fprintf(stderr, "\n");
    
    // close file
    fclose(fp);
    
    printf("Total seqs processed: %lu\tLength: %zu\tChunk: %zu\n", totalN, LEN, CHUNK_SZ);
    print_pos_stats(d_pos_stats, LEN, totalN);
    cudaDeviceSynchronize();    
    
    // cleanuup
    // device
    cudaFree(d_str);
    cudaFree(d_pos_stats);
    
    // host
    free(h_str);
    
    return 0;
}