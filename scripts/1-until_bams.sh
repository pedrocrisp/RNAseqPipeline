#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/../"

source "$basedir/common.sh"

wsdir="/home/kevin/ws"
refdir="${wsdir}/refseqs"

sample=$1

# setup

if [ ! -d reads ]
then
	echo "No ./reads/ directory"
	exit -1
fi

if [ ! -d reads/${sample} ]
then
	echo "Error: sample '${sample}' does not exist"
	exit -1
fi

# QC
echo "Initial FastQC"
mkdir -p qc/before/${sample}
time bash ${basedir}/01-qc/fastqc.sh -i reads/${sample} -o qc/before/${sample} -a ""
#check_return

echo "Run scythe"
mkdir -p qcd/${sample}/
time bash ${basedir}/01-qc/scythe.sh -i reads/${sample} -o qcd/${sample} -a "-p 0.1 -a ${refdir}/adaptors.fa"
#check_return

echo "Run FastQC after all QC steps"
mkdir -p qc/after/${sample}
time bash ${basedir}/01-qc/fastqc.sh -i qcd/${sample} -o qc/after/${sample} -a "-q"
#check_return


# align
echo "Align with subread"
mkdir -p align/${sample}
time bash ${basedir}/02-align/subread.sh -i qcd/${sample} -o align/${sample} -a "-i ${refdir}/TAIR10_gen/TAIR10_gen"
#check_return
