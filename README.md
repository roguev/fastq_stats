# fastq_stats
This repo contains a small programs in C for computing basic sequence statistics from FASTQ files.

fastq_stats.cc is a serial implementation in C.

cu_fastq_stats.cu and cu_fastq_stats_nmap.cu are parallel implementations in CUDA C differring by the way the FASTQ file is read. 

For small files (up to about 10M reads) the CUDA implementations are about twice as fast as the serial one and are able to process ca. 1M reads / sec using a GTX 580 graphics card. For larger files the performance is pretty much the same because of a HDD I/O bottleneck.
