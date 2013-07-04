#!/bin/bash

baseDistDir="$(dirname $(readlink -f $0))"
isSubDir="$(basename $baseDistDir | grep -P "^\d+-")"
if [ -n $isSubDir ]
then
	baseDistDir="${baseDistDir}/../"
fi
source "${baseDistDir}/common.sh"

getDefaultOptions $@
