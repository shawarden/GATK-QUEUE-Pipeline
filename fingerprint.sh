#!/bin/bash

# This scrip will parse a vcf file and replace the individual ID
# 
# Once generated it can be used to compare any other samples from this platform.
#
# @params
#  i  Input VCF file
#  f  Generic Identity (defaults to "FingerPrint"
#  o  Output VCF file
#

echo "Generate generic fingerprint VCF that can be compared to any other."
echo ""
echo "INFO: Commandline:      ${0} ${@}"
echo ""


usage() {
cat << EOF
This scrip will parse a vcf file and replace the individual ID

 Once generated it can be used to compare any other samples from this platform.

usage: $0 options

Options:
  -i Source VCF file
  -f Generic Identity
  -o Ouput VCF file
EOF
}

INPUT_FILE=
FINGERPRINT="Individual"
OUTPUT_FILE=

while getopts "i:f:o:" OPTION
do
	case $OPTION in
		i)
			if [[ -e ${OPTARG} ]]
			then
				INPUT_FILE=${OPTARG}
			else
				echo "FAILURE: Input file \"${OPTARG}\" does not exist."
			fi
			;;
		f)
			FINGERPRINT=${OPTARG}
			;;
		o)
			if [[ -e ${OPTARG} ]]
			then
				echo "WARNING: Output file \"${OPTARG}\" already exists and will be overwritten!"
			fi
			
			OUTPUT_FILE=${OPTARG}
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [[ -z ${INPUT_FILE} ]] || [[ -z ${OUTPUT_FILE} ]]
then
	usage
	exit 1
fi

echo Input VCF File:  ${INPUT_FILE}
echo Output VCF File: ${OUTPUT_FILE}
echo Generic Identity:${FINGERPRINT}

OUTPUT_DIR=$(dirname ${OUTPUT_FILE})
if mkdir -p ${OUTPUT_DIR}
then
	if sed '/#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT/c\#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	'${FINGERPRINT} ${INPUT_FILE} > ${OUTPUT_FILE}
	then
		echo "SUCCESS: Fingerprint file \"${OUTPUT_FILE}\" created with Generic ID \"${FINGERPRINT}\""
		exit 0
	fi
fi

echo "FAILURE: Unable to replace identifier in \"${INPUT_FILE}\""

exit 1
