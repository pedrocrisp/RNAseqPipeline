#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

bamfile="$(findBAMFiles ${input})"
if [ -z "${bamfile}" ]
then 
	samfile="$(findSAMFiles ${input})"
	# sam
	echo "featureCounts -i $samfile -o \"$output/$(basename $output).counts\" $args"
	featureCounts -i $samfile -o "$output/$(basename $output).counts" $args
else
	# bam
	echo "featureCounts -b -i $bamfile -o \"$output/$(basename $output).counts\" $args"
	featureCounts -b -i $bamfile -o "$output/$(basename $output).counts" $args
fi
