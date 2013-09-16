#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

sample=$(basename $input)
fqFiles="$(ls $input/*.f[aq]*)"
numFqFiles=$(echo $fqFiles | wc -w)

outsam="${output}/${sample}.sam"
outbam="${output}/${sample}" # no .bam, as samtools sort -f has a bug.
tmpbam="${output}/${RANDOM}.bam"

if [ ${numFqFiles} -eq 1 ]
then
	echo subjunc $args -r "$fqFiles" -o "$outsam"
	subjunc $args -r "$fqFiles" -o "$outsam"
elif [ ${numFqFiles} -eq 2 ]
then
	fq1="$(echo $fqFiles |cut -d ' ' -f 1)"
	fq2="$(echo $fqFiles |cut -d ' ' -f 2)"
	echo subjunc $args -r "${fq1}" -R "${fq2}" -o "$outsam"
	subjunc $args -r "${fq1}" -R "${fq2}" -o "$outsam"
else
	echo "ERROR: not able to align multiple fq files per pair"
	echo "fqFiles:"
	echo "${fqFiles}"
	exit 1
fi

echo "samtools view -S -u $outsam > ${tmpbam}
samtools sort -m 2G ${tmpbam} $outbam
samtools index ${outbam}.bam
rm -v ${outsam} ${tmpbam}"

samtools view -S -u $outsam > ${tmpbam}
samtools sort -m 2G ${tmpbam} $outbam
samtools index ${outbam}.bam
rm -v ${outsam} ${tmpbam}samtools view -S -u $outsam > ${tmpbam}
