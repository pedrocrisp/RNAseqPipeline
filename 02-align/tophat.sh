#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

# order of arguments to tophat.sh is important. Ref base must be given to -a, without an argument, and must be the last argument
getDefaultOptions $@

sample=$(basename $input)
fqFiles="$(findFastqFiles $input)"
numFqFiles=$(echo $fqFiles | wc -w)

outbam="${output}/${sample}" # no .bam, as samtools sort -f has a bug.
tophatbam="${output}/tophat_out/accepted_hits.bam"

if [ ${numFqFiles} -eq 1 ]
then
	echo tophat --keep-fasta-order --rg-library ${sample} -o "$output/tophat_out" $args "$fqFiles"
	tophat --keep-fasta-order --rg-library ${sample} -o "$output/tophat_out" $args "$fqFiles"
elif [ ${numFqFiles} -eq 2 ]
then
	fq1="$(echo $fqFiles |cut -d ' ' -f 1)"
	fq2="$(echo $fqFiles |cut -d ' ' -f 2)"
	echo tophat --keep-fasta-order --rg-library ${sample} -o "$output/tophat_out" $args "${fq1}" "${fq2}"
	tophat --keep-fasta-order --rg-library ${sample} -o "$output/tophat_out" $args "${fq1}" "${fq2}"
else
	echo "ERROR: not able to align multiple fq files per pair"
	echo "fqFiles:"
	echo "${fqFiles}"
	exit 1
fi


echo "ln -s $tophatbam $outbam
samtools index ${outbam}.bam"

ln -s $tophatbam $outbam
samtools index ${outbam}.bam
