/* ===================================================================
 * fastq_stats.cc
 *  computes basic sequence stats for a FASTQ dataset
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

/* reads a single FASTQ sequence record from a file
 * FASTQ record consiste of 4 lines
 *  header starting with '@'
 *  sequence string
 *  extra string starting with '+'
 *  quality string
 * The function stitches the sequence and quality strings together
 * into the output string
 * Arguments:
 *  char* to hold the result
 *  FILE* to an open file
 *  size_t maximum to length to read
*/
bool read_fastq(char* str, FILE* fp, const size_t L) {
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
                memcpy(str,line,rs*sizeof(char)); 
            else
                memcpy(str,line,L*sizeof(char));
            
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
                memcpy(str+L,line,rs*sizeof(char));
            else 
                memcpy(str+L,line,L*sizeof(char));
            
            free(line);
            return true;    // success
        
        } else
            break;
    }
    
    // something went wrong
    free(line);
    return false;
}

/* computes the sequence stats 
 * result goes in a linear array with stride 6
 * for A,T,C,G,N and quality stats
 * Arguments:
 *  unint64_t* to hold the resule
 *  char* to read from
 *  size_t length to read
 * */
void calc_stats(float* stats, const char* str, const size_t L) {
    // per-position base compositon
    for (size_t i = 0; i < L; i++) {
        switch(str[i]) {
            case 'A': case 'a':
                stats[6*i+0]++;
                break;
            case 'T': case 't':
                stats[6*i+1]++;
                break;
            case 'G': case 'g':
                stats[6*i+2]++;
                break;
            case 'C': case 'c':
                stats[6*i+3]++;
                break;
            case 'N': case 'n':
                stats[6*i+4]++;
                break;
        }
        // per-position quality
        stats[6*i+5] += str[L+i];
    }
}

/* prints the stats array
 * Arguents:
 *  unint64_t* input stats srray
 *  size_t to length
 *  unint64_t to number of sequences used for averaging
 */
void print_stats(float* stats, const size_t L, const uint64_t N) {
    const char charset[] = "ATGCNQ";
    for (int i = 0; i < 6; i++) {
        printf("@%c", charset[i]);
        for (size_t j = 0; j < L; j++) {
            // compute averages
            if (i < 5)
                printf("\t%0.4f", (float) stats[6*j + i]/N);
            else {
                if (stats[6*j + i] == 0)
                    printf("\t%0.4f", (float)0.0);
                else
                    printf("\t%0.4f", (float)(stats[6*j + i] - 33*N)/N);
            }
        }
        printf("\n");
    }
}

/* expects 3 command line arguments
 *  file to read from
 *  maximum length
 *  maximum number of sequences to process - useful to just get a rough estimate
 */ 
int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "%s: Insufficient arguments\n", argv[0]);
        exit(1);
    }
    
    char* fname = argv[1];              // input filename
    size_t L = (size_t)atoi(argv[2]);   // desired sequence lenght
    
    // max esired sequences, from command line
    uint64_t maxN = 0;
    if (argc == 4) { maxN = (uint64_t)atoi(argv[3]); }  // set if specified
    
    uint64_t N = 0;     // sequence counter
 
    // open file for reading
    FILE* fp = fopen(fname, "r");
    if (fp == NULL) {
        err(1, "fopen: %s", fname);
        exit(2);
    }
    
    // allocate space for data on the heap
    // free later
    char* str = (char*)calloc(2*L, sizeof(char));
    float* stats = (float*)calloc(6*L, sizeof(float));
    
    // read data from file
    while(1) {
        // reset string memory
        memset(str, 0, 2*L*sizeof(char));
        
        // check if reached maxN. ignore if maxN == 0
        if ((maxN != 0) && (N == maxN))
            break;
        
        // read a sequence and calc stats
        if (read_fastq(str, fp, L)) {
            calc_stats(stats,str,L);
            N++;
        } else
            break;
    }
    fclose(fp);     // close file
    
    // output
    printf("Seqs: %lu\n", N);
    print_stats(stats, L, N);
 
    // cleanup
    free(stats);
    free(str);
    return 0;
}
