#!/bin/bash

# shellcheck source=scripts/utils.sh
source "${LG3_HOME:?}/scripts/utils.sh"
source_lg3_conf

PROGRAM=${BASH_SOURCE[0]}
echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] BEGIN: $PROGRAM"
echo "Call: ${BASH_SOURCE[*]}"
echo "Script: $PROGRAM"
echo "Arguments: $*"

FAILED=false

### Configuration
LG3_HOME=${LG3_HOME:?}
LG3_OUTPUT_ROOT=${LG3_OUTPUT_ROOT:-output}
LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT:?}
LG3_DEBUG=${LG3_DEBUG:-true}

### Debug
if [[ $LG3_DEBUG ]]; then
  echo "LG3_HOME=$LG3_HOME"
  echo "LG3_OUTPUT_ROOT=$LG3_OUTPUT_ROOT"
  echo "LG3_SCRATCH_ROOT=$LG3_SCRATCH_ROOT"
  echo "PWD=$PWD"
  echo "USER=$USER"
fi


RED='\033[0;31m'
GRN='\033[0;32m'
YEL='\033[0;33m'
NOC='\033[0m'

if [ $# -lt 2 ]; then
    error "please specify project and patient!"
fi

PROJECT=$1
PATIENT=$2
CONV=${CONV:-patient_ID_conversions.tsv}
WORKDIR=${LG3_OUTPUT_ROOT}/${PROJECT:?}/mutations/${PATIENT}_mutect

echo -e "Checking MuTect output for ${YEL}${PATIENT}${NOC}, project ${PROJECT}"
echo "conversion ${CONV}"

## Pull out patient specific conversion info
grep -P "\\t${PATIENT}\\t" "${CONV}" | tr -d '\r' > "${PATIENT}.temp.conversions.txt"

## Get normal ID
while IFS=$'\t' read -r ID _ _ SAMP
do
        if [ "${SAMP:0:6}" = "Normal" ]; then
                normid=${ID}
                break
        fi
done < "${PATIENT}.temp.conversions.txt"

echo "- normid='${normid:?}'"

## Cycle through tumors and submit MUTECT jobs
while IFS=$'\t' read -r ID _ _ SAMP
do
        if [ "${SAMP:0:6}" = "Normal" ]; then
                continue
        elif [ "${SAMP:0:2}" = "ML" ]; then
                samp_label="ML"
        elif [ "${SAMP:0:3}" = "GBM" ]; then
                samp_label="GBM"
        elif [ "${SAMP:0:3}" = "Pri" ]; then
                samp_label="TUM"
        elif [ "${SAMP:0:3}" = "Tum" ]; then
                samp_label="TUM"
        elif [ "${SAMP:0:11}" = "Recurrence1" ]; then
                samp_label="REC1"
        elif [ "${SAMP:0:11}" = "Recurrence2" ]; then
                samp_label="REC2"
        elif [ "${SAMP:0:11}" = "Recurrence3" ]; then
                samp_label="REC3"
        elif [ "${SAMP:0:5}" == "tumor" ]; then
                samp_label="unkTUM"
        else
                samp_label="TUM"
        fi

        echo "- ID='${ID:?}'"
        echo "- samp_label='${samp_label:?}'"
	
        ## Expected output:
        OUT=$WORKDIR/${PATIENT}.NOR-${normid}__${samp_label}-${ID}.annotated.mutations
OK="$GRN OK$NOC"
ERR="$RED missing$NOC"
        if [ -s "$OUT" ]; then
                echo -e "$ID $OK"
        else
                echo -e "$ID $ERR"
                FAILED=true
        fi

done < "${PATIENT}.temp.conversions.txt"

## Delete PATIENT specific conversion file
rm "${PATIENT}.temp.conversions.txt"

${FAILED} && error "Script failed"

echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] END: $PROGRAM"

