#1/bin/bash

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
	find $1 -name \*.fq\* -or -name \*.fastq\*
}

alias check_return="if [ $? -ne 0] ; then exit $?; fi"
