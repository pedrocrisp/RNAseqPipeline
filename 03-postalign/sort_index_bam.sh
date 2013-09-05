#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

sams=$(find $input -type f -name *.sam -not -name *.sorted.sam)

for sam in $sams
do
	sample="$(basename $sam)"
	outbam="${output}/${sample}.sorted.bam"
	samtools view -Su $sam | samtools sort -m 2G -f - $outbam && \
		samtools index $outbam
done

bams=$(find $input -type f -name *.bam -not -name *.sorted.bam)
for bam in $bams
do
	sample="$(basename $bam)"
	outbam="${output}/${sample}.sorted.bam"
	samtools sort -m 2G -f $bam $outbam && \
		samtools index $outbam
done
