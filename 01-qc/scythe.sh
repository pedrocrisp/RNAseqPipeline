#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"
getDefaultOptions $@

for fq in $input/*.f[aq]*
do
	echo "$fq"
	sample=$(basename $output)
	outputFile="$output/${sample}.trimmed.fastq"
	echo scythe $args -o $outputFile $fq
	scythe $args -o $outputFile $fq
	gzip $outputFile
done
