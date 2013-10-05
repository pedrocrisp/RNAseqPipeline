#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

sample=$(basename $input)
fqFiles="$(findFastqFiles $input)"
numFqFiles=$(echo $fqFiles | wc -w)

outsam="${output}/${sample}.sam"
outbam="${output}/${sample}.bam"

if [ ${numFqFiles} -eq 1 ]
then
	echo subread-align $args -r "$fqFiles" -o "$outsam"
	subread-align $args -r "$fqFiles" -o "$outsam"
elif [ ${numFqFiles} -eq 2 ]
then
	fq1="$(echo $fqFiles |cut -d ' ' -f 1)"
	fq2="$(echo $fqFiles |cut -d ' ' -f 2)"
	echo subread-align $args -r "${fq1}" -R "${fq2}" -o "$outsam"
	subread-align $args -r "${fq1}" -R "${fq2}" -o "$outsam"
else
	echo "ERROR: not able to align multiple fq files per pair"
	echo "fqFiles:"
	echo "${fqFiles}"
	exit 1
fi

echo "samtools view -S -u $outsam > ${outbam}
rm -v ${outsam}"

samtools view -S -u $outsam > ${outbam}
rm -v ${outsam}
