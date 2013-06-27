#!/bin/bash
# Get dir this script is stored in, so we can call others
# FOR DEBUGGING: ScriptDir="~/Documents/pete-pipeline"
ScriptDir="$( cd "$( dirname "$0" )" && pwd )"


NTHREADS=1 # Default number of threads, in case we forget to updateNThreads
function updateNThreads {
	# Get the number of cores to use. = num processeors - 1-min load average - 1.
	NTHREADS=$(($(nproc) - $(uptime |grep -oP ": \d+\." |sed -e "s/^\: //" |sed -e "s/\.$//") - 1))
	echo "Can use $NTHREADS threads"
}

# enter directory containing reads
cd /media/work/pete/Project_SN700819R_0100_PCrisp_RSB_Arabidopsis_mRNA/

### DO QC
mkdir QC

# run fastqc
updateNThreads
fastqc Sample_*/*.fastq.gz -o QC -t $NTHREADS


# and view the results (remember Ctrl-tab changes tabs)
for qc_res in QC/BJP_*_fastqc/fastqc_report.html
do
	firefox "file://$(readlink -f $qc_res)" &
done

# get full path to base sample directories on
workingdir=$(pwd)

# Do trim_galore
updateNThreads
# Run TG
	# min qual of 20
	# run fastqc on trimmed file
	# min length 20
	# Compress the output file with GZIP
	# use pared-end
	# run on the *R1*.fastq.gz and *R2*.fastq.gz files

find -maxdepth 1  -type d  -name Sample_\* | parallel -j $NTHREADS "cd {}; trim_galore \
	-q 20 \
	--fastqc \
	--length 20 \
	--gzip \
	--paired \
	*.fastq.gz \
	>./trim_galore.log 2>&1"
