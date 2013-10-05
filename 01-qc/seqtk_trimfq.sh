#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"
getDefaultOptions $@

for fq in $(findFastqFiles $input/)
do
	fqname="$(basename $fq)"
	sample=$(basename $output)
	outputFile="$output/${fqname%%.*}.trimmed.fq"
	echo "seqtk trimfq $args $fq >$outputFile"
	seqtk trimfq $args $fq >$outputFile
done
