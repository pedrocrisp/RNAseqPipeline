#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"
getDefaultOptions $@

scytheOut="$(basename $output .gz)"
echo scythe $args -o $output $input
scythe $args -o $scytheOut $input
gzip $scytheOut
