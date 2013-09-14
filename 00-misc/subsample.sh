#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"
getDefaultOptions $@

for fq in $(find $input/ -name *.f[aq]*)
do
	fqname="$(basename $fq)"
	sample=$(basename $output)
	outputFile="$output/${fqname%%.*}.$args.fq.gz"
	echo "seqtk sample \"$fq\" $args | gzip - >\"$outputFile\""
	seqtk sample "$fq" $args | gzip - >"$outputFile"
done
