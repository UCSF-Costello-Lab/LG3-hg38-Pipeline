#!/bin/bash

# shellcheck source=scripts/utils.sh
source "${LG3_HOME:?}/scripts/utils.sh"
source_lg3_conf

PROGRAM=${BASH_SOURCE[0]}
PROG=$(basename "$PROGRAM")
echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] BEGIN: $PROGRAM"
echo "Call: ${BASH_SOURCE[*]}"
echo "Script: $PROGRAM"
echo "Arguments: $*"

### Configuration
LG3_HOME=${LG3_HOME:?}
LG3_OUTPUT_ROOT=${LG3_OUTPUT_ROOT:-output}
LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT:?}
LG3_DEBUG=${LG3_DEBUG:-true}
ncores=${SLURM_NTASKS:-1}
VERBOSITY=ERROR
CLEAN=true

TMP_DIR="${LG3_SCRATCH_ROOT}/tmp"
make_dir "${TMP_DIR}"

### Debug
if $LG3_DEBUG ; then
  echo "$PROG Settings:"
  echo "- LG3_HOME=$LG3_HOME"
  echo "- LG3_OUTPUT_ROOT=$LG3_OUTPUT_ROOT"
  echo "- LG3_SCRATCH_ROOT=$LG3_SCRATCH_ROOT"
  echo "- TMP_DIR=$TMP_DIR"
  echo "- PWD=$PWD"
  echo "- USER=$USER"
  echo "- hostname=$(hostname)"
  echo "- ncores=$ncores"
fi

#module load gatk/4.3.0.0
echo "GATK4 = ${GATK4}"
module load openjdk/1.8.0
echo "Java : "
java -version 2>&1

## Input
BAM=$1
PATIENT=$2
ILIST=$3
echo "Input:"
echo "- BAM=${BAM:?}"
echo "- PATIENT=${PATIENT:?}"
echo "- ILIST=${ILIST:?}"

## Assert existance of input files
assert_file_exists "${BAM}"
assert_file_exists "${ILIST}"

## References
echo "References:"
echo "- REF=${REF:?}"
assert_file_exists "${REF}"

## Software
echo "Software:"
echo "- GATK=${GATK4:?}"


## Assert existance of software
assert_file_executable "${JAVA}"
assert_file_exists "${GATK4}"

HaplotypeCaller() {
   echo -e "\\n\\n****** HaplotypeCaller : Call SNP and indel variants in ERC mode ******"
	bambase=$(basename "$1")
	ID=${bambase%%.*}
	echo "Processing ${ID} ..."
   { time ${GATK4} HaplotypeCaller \
      --verbosity ${VERBOSITY} \
		--tmp-dir "${TMP_DIR}" \
		-L "${ILIST}" \
		-ip "${PADDING}" \
      -R "${REF}" \
		-I "$1" \
      -O "${ID}.g.vcf" \
      -ERC GVCF; } 2>&1 || error "HaplotypeCaller FAILED"
   assert_file_exists "${ID}.g.vcf"
	echo -ne "Called ${ID}.g.vcf : "
	grep -cv '^#' "${ID}".g.vcf
   echo "****** HaplotypeCaller Completed! ******"
}

echo "BAM = ${BAM}"
HaplotypeCaller "${BAM}"

echo "[Germline] Part 1 Finished!"

echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] END: $PROGRAM"
