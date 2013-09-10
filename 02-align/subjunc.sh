#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

sample=$(basename $input)
fqFiles="$(ls $input/*.f[aq]*)"
numFqFiles=$(echo $fqFiles | wc -w)

outsam="align/${sample}/${sample}.sam"
outbam="align/${sample}/${sample}.bam"
tmpbam="align/${sample}/${RANDOM}.bam"

if [ ${numFqFiles} -eq 1 ]
then
	subread-align $args -r "$fqFiles" -o "$outsam"
elif [ ${numFqFiles} -eq 2 ]
then
	fq1="$(echo $fqFiles |cut -d ' ' -f 1)"
	fq2="$(echo $fqFiles |cut -d ' ' -f 2)"
	subjunc $args -r "${fq1}" -R "${fq2}" -o "$outsam"
else
	echo "ERROR: not able to align multiple fq files per pair"
	echo "fqFiles:"
	echo "${fqFiles}"
fi

samtools view -S -u $outsam > ${tmpbam}
samtools sort -m 2G -f ${tmpbam} $outbam
samtools index $outbam
rm -v ${outsam} ${tmpbam}