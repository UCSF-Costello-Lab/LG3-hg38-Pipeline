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
LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT:?}
LG3_DEBUG=${LG3_DEBUG:-true}

### Debug
if [ $LG3_DEBUG = true ]; then
  echo "Settings:"
  echo "- LG3_HOME=${LG3_HOME:?}"
  echo "- LG3_OUTPUT_ROOT=${LG3_OUTPUT_ROOT:?}"
  echo "- LG3_SCRATCH_ROOT=$LG3_SCRATCH_ROOT"
  echo "- PWD=$PWD"
  echo "- USER=$USER"
  echo "  ncores=${SLURM_NTASKS}"
  echo "- hostname=$(hostname)"
  echo "- TG=${TG:?}"
  echo "- CUTADAPT=${CUTADAPT:?}"
fi

assert_python ""

QTY=20
LEN=20
STRINGENCY=1

if [ $# -lt 3 ]; then
    error "Run trim_galore in paired mode on gzipped fastq files using Illumina universal adapter\\nUsage: $0 [ -quality=$QTY -length=$LEN -stringency=$STRINGENCY] 1.fastq.gz 2.fastq.gz sampleID"
fi

#### Parse optional args
while [ -n "$1" ]; do
case $1 in
    -q*=*) QTY=${1#*=};shift 1;;
    -l*=*) LEN=${1#*=};shift 1;;
    -s*=*) STRINGENCY=${1#*=};shift 1;;
    -*) error "No such option $1";;
    *)  break;;
esac
done

FQ1=$1
FQ2=$2
SAMPLE=$3

### Input
echo "Input:"
echo "FQ1=${FQ1:?}"
echo "FQ2=${FQ2:?}"
echo "LEN=${LEN:?}"
echo "QTY=${QTY:?}"
echo "STRINGENCY=${STRINGENCY:?} (ignored; always STRINGENCY=1)"

## Assert existance of input files
assert_file_exists "${FQ1}"
assert_file_exists "${FQ2}"

echo "Software:"
echo "- TG=${TG:?}"
assert_file_executable "${TG}"
echo "- CUTADAPT=${CUTADAPT:?}"

### Default: --length 20 --quality 20 --stringency 1
### --fastqc
time $TG --paired --quality "$QTY" --length "$LEN" --stringency 1 --path_to_cutadapt "$CUTADAPT" --illumina "$FQ1" "$FQ2" || error "[trim_galore] trim_galore FAILED"

## Assert existance of output files
assert_file_exists "${SAMPLE}_R1"*_val_1.fq.gz
assert_file_exists "${SAMPLE}_R2"*_val_2.fq.gz
assert_file_exists "${SAMPLE}_R1"*.fastq.gz_trimming_report.txt
assert_file_exists "${SAMPLE}_R2"*.fastq.gz_trimming_report.txt

echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] END: $PROGRAM"

