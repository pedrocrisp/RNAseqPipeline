#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/.."

source "$basedir/common.sh"
getDefaultOptions $@

fastqc $args -o "$output" $input/*.f[aq]*
