#!/bin/bash

NTHREADS=23 # Default number of threads, in case we forget to updateNThreads
function updateNThreads {
	# Get the number of cores to use. = num processeors - 1-min load average - 1.
	NTHREADS=$(($(nproc) - $(uptime |grep -oP ": \d+\." |sed -e "s/^\: //" |sed -e "s/\.$//") - 1))
	echo "Can use $NTHREADS threads"
}

# enter directory containing reads
cd /home/pete/workspace/Project_SN700819R_0100_PCrisp_RSB_Arabidopsis_mRNA/mapping

#number of threads to use
#updateNThreads

# Run samtools view in parallel
	
find -maxdepth 2  -type d  -name Sample_\* | parallel -j $NTHREADS "cd {}; \
	samtools sort -n -m 3G tophat_out/accepted_hits.bam tophat_out/accepted_hits_sorted; samtools view \
	-h \
	-o tophat_out/accepted_hits_sorted.sam \
	tophat_out/accepted_hits_sorted.bam \
	>./samtools.log 2>&1"

