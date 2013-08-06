#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/../"

source "$basedir/common.sh"


###### setup #######
wsdir="/home/kevin/ws"
refdir="${wsdir}/refseqs"

sample=$1


##### Check env ####

if [ ! -d reads ]
then
	echo "No ./reads/ directory"
	exit -1
fi

if [ -z "${sample}" ]
then
	echo "Must give sample"
fi

if [ ! -d "reads/${sample}" ]
then
	echo "Error: sample '${sample}' does not exist"
	exit -1
fi

######## QC ########
echo "Initial FastQC"
mkdir -p qc/before/${sample}
time bash ${basedir}/01-qc/fastqc.sh -i reads/${sample} -o qc/before/${sample} -a ""

echo "Run scythe"
qcstep="scythe"
mkdir -p qcd/${qcstep}/${sample}/
time bash ${basedir}/01-qc/scythe.sh -i reads/${sample} -o qcd/${qcstep}/${sample} -a "-p 0.1 -a ${refdir}/adaptors.fa"

echo "Run seqtk"
qcstep="trimfq"
mkdir -p qcd/${qcstep}/${sample}/
time bash ${basedir}/01-qc/seqtk_trimfq.sh -i reads/${sample} -o qcd/${qcstep}/${sample} -a ""

# enter any additional qc here

pushd qcd > /dev/null
ln -s ${qcstep}/${sample} ${sample}
popd >/dev/null

echo "Run FastQC after all QC steps"
mkdir -p qc/after/${sample}
time bash ${basedir}/01-qc/fastqc.sh -i qcd/${sample} -o qc/after/${sample} -a ""


###### align #######
echo "Align with subread"
mkdir -p align/${sample}
time bash ${basedir}/02-align/subread.sh -i qcd/${sample} -o align/${sample} -a "-i ${refdir}/TAIR10_gen/TAIR10_gen"

###### count #######
echo "Count with featurecounts"
mkdir -p count/${sample}
time bash ${basedir}/04-initialstats/featurecounts.sh -i align/${sample} -o count/${sample} -a "-F SAF -b -a ${refdir}/TAIR10_gen/TAIR10_GFF3_genes_transposons.tab -p -C"
