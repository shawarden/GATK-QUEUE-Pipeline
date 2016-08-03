#!/bin/bash
export SCRIPT_PATH="/home/clinical_genetics/bin/q-shell/FastQtoGVCF.q"
export RUN_NAME="FQ2GVCF"

if /home/clinical_genetics/bin/q-shell/queue ${*}
then
	exit 0
else
	exit 1
fi
