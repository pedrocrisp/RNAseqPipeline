#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"
getDefaultOptions $@

fqfiles=$(findFastqFiles "$input")
fastqc $args -o "$output" $fqfiles
