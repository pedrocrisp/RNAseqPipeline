#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

if [ ${numFqFiles} -eq 1 ]
then
	tophat $args  -o "$output/tophat_out" "$fqFiles"
elif [ ${numFqFiles} -eq 2 ]
then
	fq1="$(echo $fqFiles |cut -d ' ' -f 1)"
	fq2="$(echo $fqFiles |cut -d ' ' -f 2)"
	tophat $args -o "$output/tophat_out" "${fq1}" "${fq2}"
else
	echo "ERROR: not able to align multiple fq files per pair"
	echo "fqFiles:"
	echo "${fqFiles}"
fi

# link to keep bam position compatible w/ other aligners
ln -s "$output/tophat_out/accepted_hits.bam" "$output/$(basename $output).bam"
