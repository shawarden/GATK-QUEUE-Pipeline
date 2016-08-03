#!/bin/bash

# This shell script is used to connect library merge to individual merging.
#
# Takes a file (presumably a bam file) creates a hard link to it.
#

echo "Pipeline junction between merging to Library and merging to Individual."
echo ""
echo "INFO: Commandline:      ${0} ${@}"
echo ""

usage() {
cat << EOF
 This shell script is used to connect library merge to individual merging.

 Takes a file (presumably a bam file) creates a hard link to it.
  

usage: $0 -b <file.bam> -i <file.bai> -l <file.bam> -o <file.bai>

Options:
  -b Source BAM file
  -i Source BAI file
  -l Output BAM file
  -o Output BAI file
EOF
}

ORIGIN_BAM_FILE=
ORIGIN_BAI_FILE=
LINKED_BAM_FILE=
LINKED_BAI_FILE=

while getopts "b:i:l:o:" OPTION
do
	case $OPTION in
		b)	# BAM file. Has to exist!
			if [[ -e ${OPTARG} ]]
			then
				ORIGIN_BAM_FILE=${OPTARG}
			else
				echo "FAILURE: BAM file \"${OPTARG}\" does not exist."
			fi
			;;
		i)	# BAI file. Has to exist!
			if [[ -e ${OPTARG} ]]
			then
				ORIGIN_BAI_FILE=${OPTARG}
			else
				echo "FAILURE: BAM file \"${OPTARG}\" does not exist."
			fi
			;;
		l)	# Linked BAM file. Probably shouldn't exist but not a hugh deal if it does.
			if [[ -e ${OPTARG} ]]
			then
				echo "WARNING: Link target file \"${OPTARG}\" already exists!"
				if rm ${OPTARG}
				then
					echo "WARNING: Deleting previous \"${OPTARG}\" file!"
				else
					echo "FAILURE: Unable to delete existing \"${OPTARG}\"."
					exit 1
				fi
			fi
			
			LINKED_BAM_FILE=${OPTARG}
			;;
		o)	# Linked BAI file. Probably shouldn't exist but not a hugh deal if it does.
			if [[ -e ${OPTARG} ]]
			then
				echo "WARNING: Link target file \"${OPTARG}\" already exists!"
				if rm ${OPTARG}
				then
					echo "WARNING: Deleting previous \"${OPTARG}\" file!"
				else
					echo "FAILURE: Unable to delete existing \"${OPTARG}\"."
					exit 1
				fi
			fi
			
			LINKED_BAI_FILE=${OPTARG}
			;;
		\?)
			echo "FAILURE: Invalid options -${OPTARG}"
			usage
			exit 1
			;;
	esac
done

if [[ -z ${ORIGIN_BAM_FILE} ]] || [[ -z ${ORIGIN_BAI_FILE} ]] || [[ -z ${LINKED_BAM_FILE} ]] || [[ -z ${LINKED_BAI_FILE} ]]
then
	usage
	exit 1
fi

echo "INFORMA: Source files:     ${ORIGIN_BAM_FILE} & ${ORIGIN_BAI_FILE}"
echo "INFORMA: Linked files:     ${LINKED_BAM_FILE} & ${LINKED_BAI_FILE}"

echo ""

LINKED_DIR=$(dirname ${LINKED_BAM_FILE})
mkdir -p ${LINKED_DIR}

if ln ${ORIGIN_BAM_FILE} ${LINKED_BAM_FILE}
then
	echo "SUCCESS: ${ORIGIN_BAM_FILE} hard linked to ${LINKED_BAM_FILE}."
else
	echo "FAILURE: Unable to hard link ${ORIGIN_BAM_FILE} to ${LINKED_BAM_FILE}."
	exit 1
fi

if ln ${ORIGIN_BAI_FILE} ${LINKED_BAI_FILE}
then
	echo "SUCCESS: ${LINKED_BAI_FILE} hard linked to ${LINKED_BAI_FILE}."
else
	echo "FAILURE: Unable to hard link ${LINKED_BAI_FILE} to ${LINKED_BAI_FILE}."
	exit 1
fi

exit 0
