#!/bin/bash

Project="Project_SN700819R_0100_PCrisp_RSB_Arabidopsis_mRNA"
WSDir="/home/pete/workspace/"
RefDir="/home/pete/workspace/refseqs"

ProjectDir="${WSDir}/${Project}"

# Make mapping dir
cd "${WSDir}/${Project}"
mkdir "mapping"
cd mapping

NTHREADS=1 # Default number of threads, in case we forget to updateNThreads
function updateNThreads {
	# Get the number of cores to use. = num processeors - 1-min load average - 1.
	NTHREADS=$(($(nproc) - $(uptime |grep -oP ": \d+\." |sed -e "s/^\: //" |sed -e "s/\.$//") - 1))
	echo "Can use $NTHREADS threads"
}


# Make the Sample_* dirs in the mapping folder
for sample in ${ProjectDir}/Sample_*
do
	sample="$(basename "$sample")"
	WorkingDir="${WSDir}/${Project}/mapping/${sample}"
	mkdir ${WorkingDir} 
	cd "${WorkingDir}"

	SampleDir="${WSDir}/${Project}/${sample}"

	updateNThreads
	tophat2 \
		--max-intron-length 40 \
		--max-intron-length 2000 \
		-G "${RefDir}/TAIR10_gen/TAIR10_GFF3_genes_transposons.gff" \
		--transcriptome-index "${RefDir}/TAIR10_tx_from_tophat/TAIR10_tx" \
		--solexa-quals \
		--library-type fr-unstranded \
		-p ${NTHREADS} \
		-g 2 \
		${RefDir}/TAIR10_gen/TAIR10_allchr \
		${SampleDir}/*.fq.gz | tee "./tophat${sample}.log"
done
