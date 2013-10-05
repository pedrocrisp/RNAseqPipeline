#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

sams=$(findSAMFiles $input)

for sam in $sams
do
	sample="$(basename $sam)"
	outbam="${output}/${sample}.sorted.bam"
	echo "samtools view -Su $sam | samtools sort -m 2G -f - $outbam && \
		samtools index $outbam"
	samtools view -Su $sam | samtools sort -m 2G -f - $outbam && \
		samtools index $outbam
done

bams=$(findBAMFiles $input)
for bam in $bams
do
	sample="$(basename $bam)"
	outbam="${output}/${sample}.sorted.bam"
	echo "samtools sort -m 2G -f $bam $outbam && \
		samtools index $outbam"
	samtools sort -m 2G -f $bam $outbam && \
		samtools index $outbam
done
