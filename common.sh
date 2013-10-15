#!/bin/bash

function getDefaultOptions () {
	while getopts :i:o:a: flag
	do
	    case $flag in
		i)
		    input="$OPTARG"
		    ;;
		o)
		    output="$OPTARG"
		    ;;
		a)
		    args=`shift $(( ${OPTIND} - 2 )); echo "${*}"`
		    break
		    ;;
	    esac
	done

	echo "Input is $input" >/dev/null
	echo "Output is $output" >/dev/null
	echo "Args are $args" >/dev/null

	if [[ -z "$input"  || ! -r "$input" ]]
	then
		echo "input was |$input|"
		echo "Must give input and it must exist"
		exit 1
	fi

	if [[ -z "$output" ]]
	then
		echo "output was |$output|"
		echo "Must give output"
		exit 1
	fi
}

function findFastqFiles () {
	find $1 -name \*.fq\* -print0 -or -name \*.fastq\* -print0 |sed -e 's/ /\\ /g' -e 's/\x0/ /g'
}

function findBAMFiles () {
	find $1 -maxdepth 1 -name \*.bam -print0 |sed -e 's/ /\\ /g' -e 's/\x0/ /g'
}

function findSAMFiles () {
	find $1 -maxdepth 1 -name \*.sam -print0 |sed -e 's/ /\\ /g' -e 's/\x0/ /g'
}

alias check_return="if [ $? -ne 0] ; then exit $?; fi"
alias timestamp='date +%Y%m%d-%H%M%S'
