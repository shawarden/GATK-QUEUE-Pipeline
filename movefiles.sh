#!/bin/bash

IN_FILE=${1}
IN_PATH=$(dirname ${IN_FILE})

OUT_FILE=${2}
OUT_PATH=$(dirname ${OUT_FILE})

echo "INFORMA Input file:  ${IN_FILE}"
echo "INFORMA Input path:  ${IN_PATH}"

echo "INFORMA Output File: ${OUT_FILE}"
echo "INFORMA Output path: ${OUT_PATH}"

if [ -d "${IN_PATH}" ]
then
	echo "INFORMA ${IN_PATH} exists."
	if [ -e "${IN_FILE}" ]
	then
		echo "INFORMA ${IN_FILE} exists within ${IN_PATH}"
		if mkdir -p ${OUT_PATH}
		then
			echo "SUCCESS Created ${OUT_PATH} destination folder."
			find ${IN_PATH} -type f -iname "*.out" -exec cp {} ${OUT_PATH} \;
			find ${IN_PATH} -type f -iname "*.g.vcf.g*" -exec cp {} ${OUT_PATH} \;
			touch ${OUT_FILE}
			exit 0
		else
			echo "FAILURE Unable to create ${OUT_PATH}"
		fi
	else
		echo "FAILURE ${IN_FILE} does not exist."
	fi
else
	echo "FAILURE ${IN_PATH} does not exist."
fi

exit 1