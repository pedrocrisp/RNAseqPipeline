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
		    args="$OPTARG"
		    ;;
	    esac
	done
	echo "Input is $input"
	echo "Output is $output"
	echo "Args are $args"

	if [[ -z "$input"  || ! -r "$input" ]]
	then
		echo "Must give input file, and it must exist"
		exit 1
	fi

	if [[ -z "$output"  || ! -r "$output" ]]
	then
		echo "Must give output file, and it must exist"
		exit 1
	fi
}
