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
	outputFile="$output/${fqname%%.*}.noadapt.fq.gz"
	echo "scythe $args $fq >$outputFile"
	scythe $args $fq >$outputFile
done
