#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/../"

source "$basedir/common.sh"


######### Setup ################
keyfile=$1

# kefile format: (tab seperated)
#Ordinal	Sample	<factor1_name> [<factor2_name>]


##### Check env ####

if [ ! -d count ]
then
	echo "No ./count directory"
	exit -1
fi


########## Run #################

cat $0
## enter steps ##


# step 2: differential expression
mkdir -p ./de
script="${basedir}/edgeR_glm_multifactor.R"
cat $script
R -f ${script} --args $keyfile ${scriptdir}/edgeR_kmhons.R >./log/de.log
