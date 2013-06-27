#!/bin/bash

RefDir="/home/pete/workspace/refseqs"

NTHREADS=1 # Default number of threads, in case we forget to updateNThreads
function updateNThreads {
	# Get the number of cores to use. = num processeors - 1-min load average - 1.
	NTHREADS=$(($(nproc) - $(uptime |grep -oP ": \d+\." |sed -e "s/^\: //" |sed -e "s/\.$//") - 1))
	echo "Can use $NTHREADS threads"
}

# enter directory containing reads
cd /home/pete/workspace/Project_SN700819R_0100_PCrisp_RSB_Arabidopsis_mRNA/mapping

#number of threads to use
updateNThreads

# Summarise read counts per gene with htseq; in parallel
# not strand specific
# 
# mode is to keep partially multimapping reads if they have a portion that unambigiously aligns (multimapped called ambugious)
##TAIR10 parent features are 3 types - gene,transposable_element_gene,pseudogene
	
find -maxdepth 2  -type d  -name Sample_\* | parallel -j $NTHREADS "cd {}; htseq-count \
	--stranded no \
	-i ID \
	-t transposable_element_gene \
	-m intersection-nonempty \
	tophat_out/accepted_hits_sorted.sam \
	${RefDir}/TAIR10_gen/TAIR10_GFF3_genes_transposons.gff \
	> tophat_out/accepted_hits.transposable_element_gene.counts \
	2>./htseq-count.errlog"


