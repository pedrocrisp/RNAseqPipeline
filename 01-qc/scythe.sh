#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"
getDefaultOptions $@

for fq in $(find $input/ -name *.f[aq]*)
do
	sample=$(basename $output)
	outputFile="$output/${sample}.noadapt.fq.gz"
	scythe $args $fq >$outputFile
done
