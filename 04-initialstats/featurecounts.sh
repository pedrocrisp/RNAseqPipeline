#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

bamfile="$(find ${input} -name \*.bam)"
featureCounts -i "$bamfile" -o "$output/$(basename $bamfile .bam).counts" $args
