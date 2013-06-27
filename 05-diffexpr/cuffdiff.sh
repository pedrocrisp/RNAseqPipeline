#!/bin/bash

Project="Project_SN700819R_0100_PCrisp_RSB_Arabidopsis_mRNA"
WSDir="/home/pete/workspace"
RefDir="/home/pete/workspace/refseqs"
ProjectDir="${WSDir}/${Project}"
Step="03-difexp"

TmpDir="/tmp/difexp_$RANDOM"
while [ -d "$TmpDir" ]
do
	TmpDir="/tmp/difexp_$RANDOM"
done
mkdir "$TmpDir" 2>/dev/null

NTHREADS=1 # Default number of threads, in case we forget to updateNThreads
function updateNThreads {
	# Get the number of cores to use. = num processeors - 1-min load average - 1.
	NTHREADS=$(($(nproc) - $(uptime |grep -oP ": \d+\." |sed -e "s/^\: //" |sed -e "s/\.$//") - 1))
	echo "Can use $NTHREADS threads"
}

	
# Make step dir
cd "${WSDir}/${Project}"
mkdir "${Step}" 2>/dev/null
cd "${Step}"

# Programs:

function edger () {
	echo "Using edgeR"
}

function cuff () {
	echo "Using cufflinks and cuffdiff"
	WorkingDir="${WSDir}/${Project}/${Step}/cuff"
	mkdir ${WorkingDir} 2>/dev/null
	cd "${WorkingDir}"
	SampleBams=""
	SampleLabels=""

	# Run cufflinks on each sample
	CufflinksDir="cufflinks"
	mkdir "${CufflinksDir}" 2>/dev/null
	cd "${CufflinksDir}"

	for aligned_sample_dir in ${ProjectDir}/mapping/Sample_*
	do
		sample="$(basename "$aligned_sample_dir")"
		mkdir "${sample}" 2>/dev/null #within cufflinks dir
		if [ -z "$SampleLabels" ]; # -z asks if string is Zero length (empty)
		then
			# If empty, dont add extra comma
			SampleLabels="$sample"
		else
			SampleLabels="$SampleLabels,$sample" # if non-empty, just append w/ comma
		fi

		SampleBam="${aligned_sample_dir}/tophat_out/accepted_hits.bam"
		# Dont need to bother with the emptyness test, as spaces are ignored on the command line
		SampleBams="$SampleBams $SampleBam"

		echo cufflinks \
			-o "${sample}" \
			-G "$RefDir/TAIR10_gen/TAIR10_GFF3_genes_transposons.gff" \
			-b "$RefDir/TAIR10_gen/TAIR10_allchr.fa" \
			-u \
			--library-type fr-unstranded \
			"${SampleBam}" \
			>>"${TmpDir}/cufflinks.cmd"
	done

	# do cufflinks in paralell
	updateNThreads
#	cat "${TmpDir}/cufflinks.cmd" | parallel --files bash -c {}

	# Cuffdiff
	cd "${WorkingDir}"
	CuffDiffDir="cuffdiff"
	mkdir ${CuffDiffDir} 2>/dev/null
	updateNThreads
	echo cuffdiff \
		-p $NTHREADS \
		-o ${CuffDiffDir} \
		-L "$(for s in Sample_BJP_317_*; do echo $s | head -c $((`echo;  done | sort -nt_ -k4,4 |uniq |tr '\n' ',' |sed -e 's/,$/\n/')"\
		-b "$RefDir/TAIR10_gen/TAIR10_allchr.fa" \
		-u \
		-T \
		--library-type fr-unstranded \
		"$RefDir/TAIR10_gen/TAIR10_GFF3_genes_transposons.gff" \
		"${SampleBams}"
}
# acutally run it
cuff
edger
rm -rv "${TmpDir}"
