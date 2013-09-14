#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"
getDefaultOptions $@

fqfiles=$(find $input/ -name *.f[aq]* -print0 |sed -e 's/\0/ /g')

echo fastqc $args -o "$output" $fqfiles
fastqc $args -o "$output" $fqfiles
