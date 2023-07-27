#!/bin/bash

# shellcheck source=scripts/utils.sh
source "${LG3_HOME:?}/scripts/utils.sh"
source_lg3_conf

PROGRAM=${BASH_SOURCE[0]}
echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] BEGIN: $PROGRAM"
echo "Call: ${BASH_SOURCE[*]}"
echo "Script: $PROGRAM"
echo "Arguments: $*"

### Configuration
LG3_HOME=${LG3_HOME:?}
LG3_OUTPUT_ROOT=${LG3_OUTPUT_ROOT:-output}
LG3_INPUT_ROOT=${LG3_INPUT_ROOT:-${LG3_OUTPUT_ROOT}}
PROJECT=${PROJECT:?}
LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT:?}
LG3_DEBUG=${LG3_DEBUG:-true}

### Debug
if [[ $LG3_DEBUG ]]; then
  echo "Settings:"
  echo "- LG3_HOME=$LG3_HOME"
  echo "- LG3_INPUT_ROOT=$LG3_INPUT_ROOT"
  echo "- LG3_OUTPUT_ROOT=$LG3_OUTPUT_ROOT"
  echo "- LG3_SCRATCH_ROOT=$LG3_SCRATCH_ROOT"
  echo "- PWD=$PWD"
  echo "- USER=$USER"
  echo "- hostname=$(hostname)"
fi

echo "Software:"
echo "- PYTHON=${PYTHON:?}"
assert_python "$PYTHON"

PYTHON_PINDEL_SETUP=${LG3_HOME}/scripts/pindel_setup.py
assert_file_exists "${PYTHON_PINDEL_SETUP}"

PATIENT=$1
PROGECT=$2
CONV=$3
echo "Input:"
echo "- PATIENT=${PATIENT:?}"
echo "- PROGECT=${PROGECT:?}"
echo "- CONV=${CONV:?}"

assert_patient_name "${PATIENT}"
assert_file_exists "${CONV}"

$PYTHON "${PYTHON_PINDEL_SETUP}" "${PATIENT}" "${PROGECT}" "${CONV}"

echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] END: $PROGRAM"
