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
	echo "featureCounts -o \"$output/$(basename $output).counts\" $args $samfile" 
	featureCounts -o "$output/$(basename $output).counts" $args $samfile 
else
	# bam
	echo "featureCounts -o \"$output/$(basename $output).counts\" $args $bamfile"
	featureCounts -o "$output/$(basename $output).counts" $args $bamfile
fi
