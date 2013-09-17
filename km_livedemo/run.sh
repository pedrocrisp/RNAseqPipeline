#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/../"

source "$basedir/common.sh"
alias timestamp='date +%Y%m%d-%H%M%S'


######### Setup ################
keyfile=$1

# kefile format: (tab seperated)
#Ordinal	Sample	<factor1_name> [<factor2_name>]


########## Run #################
# sort keyfile. -n make the header line come at the start, if it starts with a letter
sort -o $keyfile -k1n $keyfile

function getSamples() {
	grep -iv Ordinal < $keyfile | cut -f 2
}

echo "Samples are:"
echo "$(getSamples)"

cat $0
## enter steps ##

# step 1: from raw reads until counts
mkdir ./log/1-until_counts/
cat ${scriptdir}/1-until_counts.sh
getSamples |parallel "bash ${scriptdir}/1-until_counts.sh {} >./log/1-until_counts/{}.`timestamp`.log 2>&1"

# step 2: differential expression
mkdir ./de
script="${basedir}/05-diffexpr/edgeR_exact_test_pairwise.R"
cat $script
R -f $script --args $keyfile >./log/de.`timestamp`.log
