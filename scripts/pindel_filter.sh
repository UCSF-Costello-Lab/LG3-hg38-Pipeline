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
PROJECT=${PROJECT:?}
LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT:?}
LG3_DEBUG=${LG3_DEBUG:-true}

### Debug
if [[ $LG3_DEBUG ]]; then
  echo "Settings:"
  echo "- LG3_HOME=$LG3_HOME"
  echo "- LG3_OUTPUT_ROOT=$LG3_OUTPUT_ROOT"
  echo "- LG3_SCRATCH_ROOT=$LG3_SCRATCH_ROOT"
  echo "- PWD=$PWD"
  echo "- USER=$USER"
fi

#
##
### PINDEL
###
### /path/to/pindel_filter.sh <pindel_output_all>
##
#

echo "Software:"
echo "- PYTHON=${PYTHON:?}"
assert_python "$PYTHON"

PYTHON_PINDEL_FILTER=${LG3_HOME}/scripts/pindel_filter.py
assert_file_executable "${SAMTOOLS_OLD}"
echo "SAMTOOLS_OLD ${SAMTOOLS_OLD} "
assert_file_executable "${BEDTOOLS}"
assert_file_exists "${PYTHON_PINDEL_FILTER}"

datafile=$1
#proj=$2
interval=$3
echo "Input:"
echo "- datafile=${datafile:?}"
echo "- interval=${interval:?}"

assert_file_exists "${datafile}"

### filter indels
$PYTHON "${PYTHON_PINDEL_FILTER}" "${datafile}"
assert_file_exists "${datafile}.filter"
echo -n "${datafile}.filter "
wc -l "${datafile}.filter"
head -n 10 "${datafile}.filter"

### intersect with target sequence
"${BEDTOOLS}" intersect -a "${datafile}.filter" -b "${interval}" -wa > "${datafile}.filter.intersect"
assert_file_exists "${datafile}.filter.intersect"
echo -n "${datafile}.filter.intersect "
wc -l "${datafile}.filter.intersect"
head -n 10 "${datafile}.filter.intersect"

echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] END: $PROGRAM"
