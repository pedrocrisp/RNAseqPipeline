#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/../../"

source "$basedir/common.sh"

timestamp=$(date +%Y%m%d-%H%M%S)
usage="USAGE:
run.sh <keyfile>"


######### Setup ################
keyfile=$1

# kefile format: (tab seperated)
#Ordinal	Sample	<factor1_name> [<factor2_name>]

if [ ! -r $keyfile ]
then
	echo "Must provide kefile"
	echo $usage
	exit -1
fi

if [ ! -d count ]
then
	echo "No ./count directory"
	exit -1
fi


########## Run #################

cat $0

## enter steps ##

mkdir -p de

# differential expression
script="${basedir}/05-diffexpr/edgeR_pairwise.R"
logdir="./log/de.${timestamp}"
mkdir -p $logdir
cat $script
R -f $script --args $keyfile >${logdir}/edgeR_pairwise.log 2>&1
