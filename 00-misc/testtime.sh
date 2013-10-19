#!/bin/bash

#for i in `seq 1 5`
#do
#	time $@
#done

seq=$(seq 1 5)

for i in $seq
do
	echo $i aln_subread
	time bash ~/prog/bio/honsPipeline/pipelines/aln_subread/run.sh testing.key >log/${i}.aln_subread.runlog
done
for i in $seq
do
	echo $i aln_subread_htseq
	time bash ~/prog/bio/honsPipeline/pipelines/aln_subread_htseq/run.sh testing.key >log/${i}.aln_subread_htseq.runlog
done

for i in $seq
do
	echo $i aln_tophat
	time bash ~/prog/bio/honsPipeline/pipelines/aln_tophat/run.sh testing.key >log/${i}.aln_tophat.runlog
done

for i in $seq
do
	echo $i aln_tophat_htseq
	time bash ~/prog/bio/honsPipeline/pipelines/aln_tophat_htseq/run.sh testing.key >log/${i}.aln_tophat_htseq.runlog
done

