#!/bin/bash

INPUT=${1}
OUTPUT=${2}

# Skip identical file.
[ ${INPUT} -ef ${OUTPUT} ] && exit 0

INPUT_DIR=$(dirname ${INPUT})
OUTPUT_DIR=$(dirname ${OUTPUT})

INPUT_MOUNT=$(df -P -T ${INPUT_DIR} | tail -n +2 | awk '{print $2}')
OUTPUT_MOUNT=$(df -P -T ${OUTPUT_DIR} | tail -n +2 | awk '{print $2}')

INPUT_DEVICE=$(stat -c "%d" ${INPUT_DIR})
OUTPUT_DEVICE=$(stat -c "%d" ${OUTPUT_DIR})

if [ ${INPUT_MOUNT} == "cifs" ]
then
	# Source is remote file.
	if [ ${OUTPUT_MOUNT} != "cifs" ]
	then
		# destination is local file.
		echo Copy remote input file to local output file folder.
		cp ${INPUT} ${OUTPUT_DIR}
	fi
else
	# Source is local file.
	if [ ${OUTPUT_MOUNT} == "cifs" ]
	then
		# distination is remote file.
		echo Copy local input file to remote output file folder.
	fi
fi
