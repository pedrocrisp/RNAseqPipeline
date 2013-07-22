#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"

getDefaultOptions $@

tophat --solexa-quals --library-type fr-unstranded $args -o "$output/tophat_out" "$input"
