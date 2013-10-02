#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/../../"

source "$basedir/common.sh"

alias timestamp='date +%Y%m%d-%H%M%S'
alias usage="echo 'run.sh <keyfile> <r_config_file>'"


######### Setup ################
keyfile=$1
rconffile=$2

# kefile format: (tab seperated)
#Ordinal	Sample	<factor1_name> [<factor2_name>]

if [ ! -d count ]
then
	echo "No ./count directory"
	usage
	exit -1
fi

if [ ! -r $keyfile ]
then
	echo "Must provide kefile"
	usage
	exit -1
fi

if [ ! -r $rconffile ]
then
	echo "Must provide R configuration file"
	usage
	exit -1
fi

########## Run #################

cat $0
## enter steps ##


# step 2: differential expression
mkdir -p ./de
script="${basedir}/edgeR_glm_multifactor.R"
cat $script
R -f ${script} --args $keyfile $rconffile >./log/de.`timestamp`.log
