#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@


bamfile="$(findBAMFiles ${input})"
gffFile="$(echo $args| sed -e 's/ /\n/g' | grep -A 1 -- '-a' | grep -v -- '-a')"
args=$(echo $args |sed -e 's/ /\n/g' |grep -v -- -a |grep -v -- $gffFile  | tr '\n' ' ')

echo "samtools sort -on $bamfile - | samtools view -h - |htseq-count $args - $gffFile > \"$output/$(basename $output).counts\""
samtools sort -on $bamfile - | samtools view -h - |htseq-count $args - $gffFile > "$output/$(basename $output).counts"
