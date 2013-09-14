#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"
getDefaultOptions $@

for fq in $(find $input/ -name *.f[aq]*)
do
	sample=$(basename $output)
	outputFile="$output/${sample}.trimmed.fq"
	echo "seqtk trimfq $args $fq >$outputFile"
	seqtk trimfq $args $fq >$outputFile
done
