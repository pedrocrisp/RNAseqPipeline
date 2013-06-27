OVERVIEW
========

Bash script pipelines to complete an RNAseq differential expression analysis. Each directory contains scripts with alternative implementations of various parts of the pipeline. To run, create a script which ```source```es the file containing the implementation of each step you require, and run the function with the step name. See example.sh for an example.

Developed by Kevin Murray (with contributions from Peter Crisp), 2013.

LICENSE
=======

Available under the GNU GPL v3 or later.

DESCRIPTION
===========

This is a work in progress


Steps:
- 01: quality control of raw reads (trimming, generation of QC reports, etc)
- 02: alignment of raw reads to reference (everything from raw reads to a bam alignment file)
- 03: post-processing of alignemnts (exracting counts, filtering alignments by arbirtary criteria, etc)
- 04: initial statistical analysis (normalisation, etc if required [may be done in 05])
- 05: differential expression analysis (statstical testing of differential expression)
- 06: post-analysis
