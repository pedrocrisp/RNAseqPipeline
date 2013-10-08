#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/../../"

source "$basedir/common.sh"

timestamp=$(date +%Y%m%d-%H%M%S)
alias usage="echo 'run.sh <keyfile>'"


######### Setup ################
keyfile=$1

# kefile format: (tab seperated)
#Ordinal	Sample	<factor1_name> [<factor2_name>]

if [ ! -r $keyfile ]
then
	echo "Must provide kefile"
	usage
	exit -1
fi


########## Run #################
# sort keyfile. -n make the header line come at the start, if it starts with a letter
sort -o $keyfile -k1n $keyfile

function getSamples() {
	grep -iv Ordinal < $keyfile | cut -f 2 | grep --color=never "."
}

echo "Samples are:"
echo "$(getSamples)"

cat $0

## enter steps ##

# from raw reads until counts, using tophat
script=${scriptdir}/tophat.sh
logdir="./log/tophat.${timestamp}"
mkdir -p $logdir
cat $script
getSamples |parallel bash ${script} {} \>${logdir}/{}.log 2\>\&1
