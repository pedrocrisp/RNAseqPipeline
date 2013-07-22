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
		    ;;
	    esac
	done
	echo "Input is $input"
	echo "Output is $output"
	echo "Args are $args"

	if [[ -z "$input"  || ! -r "$input" ]]
	then
		echo "Must give input and it must exist"
		exit 1
	fi

	if [[ -z "$output" ]]
	then
		echo "Must give output"
		exit 1
	fi
}
