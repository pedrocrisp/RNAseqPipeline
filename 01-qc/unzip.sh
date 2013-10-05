#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"
getDefaultOptions $@

for fq in $(findFastqFiles $input)
do
	fqname="$(basename $fq)"
	outputFile="$output/${fqname%%.*}.noadapt.fq.gz"
	echo "gunzip -c $args $fq >$outputFile"
	gunzip -c $args $fq >$outputFile
done
