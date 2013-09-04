#!/bin/bash

# source common function script
scriptdir="$(dirname $(readlink -f $0))"
basedir="$scriptdir/../"

source "$basedir/common.sh"


###### setup #######
wsdir="~/workspace"
refdir="${wsdir}/refseqs"

sample=$1

# everything's already mapped, in `./mapping.report.all/<sample>`

###### count #######
echo "Count with featurecounts"
mkdir -p count/sense/${sample}
time bash ${basedir}/04-initialstats/featurecounts.sh -i mapping.report.all/${sample} -o count/sense/${sample} -a "-F SAF -b -a ${refdir}/TAIR10_gen/TAIR10_GFF3_genes.saf -s 1"
mkdir -p count/antisense/${sample}
time bash ${basedir}/04-initialstats/featurecounts.sh -i mapping.report.all/${sample} -o count/antisense/${sample} -a "-F SAF -b -a ${refdir}/TAIR10_gen/TAIR10_GFF3_genes.saf -s 1"
