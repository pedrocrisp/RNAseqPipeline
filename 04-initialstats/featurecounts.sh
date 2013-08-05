#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

bamfile="$(find ${input} -name \*.[bs]am)"
if [ -n "$(echo $bamfile | grep sam)" ]
then
	# sam
	echo "featureCounts -i \"$bamfile\" -o \"$output/$(basename $bamfile .sam).counts\" $args"
	featureCounts -i "$bamfile" -o "$output/$(basename $bamfile .sam).counts" $args
elif [ -n "$(echo $bamfile | grep bam)" ]
then
	# bam
	echo "featureCounts -b -i \"$bamfile\" -o \"$output/$(basename $bamfile .bam).counts\" $args"
	featureCounts -b -i "$bamfile" -o "$output/$(basename $bamfile .bam).counts" $args
fi
