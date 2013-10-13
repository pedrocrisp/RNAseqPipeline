#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/../../"

source "$basedir/common.sh"

timestamp=$(date +%Y%m%d-%H%M%S)
usage="USAGE:
run.sh <keyfile> <r_config_file>'"


######### Setup ################
keyfile=$1
rconffile=$2

# kefile format: (tab seperated)
#Ordinal	Sample	<factor1_name> [<factor2_name>]

if [ ! -d count ]
then
	echo "No ./count directory"
	echo $usage
	exit -1
fi

if [ ! -r $keyfile ]
then
	echo "Must provide kefile"
	echo $usage
	exit -1
fi

if [ ! -r $rconffile ]
then
	echo "Must provide R configuration file"
	echo $usage
	exit -1
fi

########## Run #################

cat $0
## enter steps ##

mkdir -p ./de

# differential expression
script="${basedir}/edgeR_glm.R"
logdir="./log/de.${timestamp}"
mkdir -p $logdir
cat $script
R -f ${script} --args $keyfile $rconffile >${logdir}/edgeR_glm.log
