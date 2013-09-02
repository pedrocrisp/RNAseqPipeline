#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

sample="$(basename $input)"
bambase="align/${sample}/${sample}"
inbam="${bambase}.bam"
outbam="${bambase}.sorted.bam"

samtools sort -f $inbam $outbam
samtools index $outbam
