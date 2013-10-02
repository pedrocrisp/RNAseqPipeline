#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/../../"

source "$basedir/common.sh"

alias timestamp='date +%Y%m%d-%H%M%S'
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

if [ ! -d count ]
then
	echo "No ./count directory"
	exit -1
fi


########## Run #################

cat $0

## enter steps ##

# step 2: differential expression
mkdir ./de
script="${basedir}/05-diffexpr/edgeR_exact_test_pairwise.R"
cat $script
R -f $script --args $keyfile >./log/de.`timestamp`.log
