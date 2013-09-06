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
	seqtk trimfq $args $fq >$outputFile
done
