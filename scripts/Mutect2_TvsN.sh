#!/bin/bash

# shellcheck source=scripts/utils.sh
source "${LG3_HOME:?}/scripts/utils.sh"
source_lg3_conf
CLEAN=true
VERBOSITY=ERROR ## WARNING INFO DEBUG
echo "- Verbosity=${VERBOSITY}"

PROGRAM=${BASH_SOURCE[0]}
echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] BEGIN: $PROGRAM"
echo "Call: ${BASH_SOURCE[*]}"
echo "Script: $PROGRAM"
echo "Arguments: $*"

#module load gatk/4.3.0.0
echo "GATK4 = ${GATK4}"
module load openjdk/1.8.0
echo "Java : "
java -version

### Configuration
LG3_HOME=${LG3_HOME:?}
LG3_OUTPUT_ROOT=${LG3_OUTPUT_ROOT:-output}
PROJECT=${PROJECT:?}
LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT:?}
LG3_DEBUG=${LG3_DEBUG:-true}
ncores=${SLURM_NTASKS:-1}

TMP_DIR="${LG3_SCRATCH_ROOT}/tmp"
make_dir "${TMP_DIR}"

### Debug
if $LG3_DEBUG ; then
  echo "Settings:"
  echo "- LG3_HOME=${LG3_HOME}"
  echo "- LG3_OUTPUT_ROOT=${LG3_OUTPUT_ROOT}"
  echo "- LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT}"
  echo "- TMP_DIR=${TMP_DIR}"
  echo "- PWD=$PWD"
  echo "- USER=${USER}"
  echo "- hostname=$(hostname)"
  echo "- ncores=${ncores}"
fi

#module load openjdk/1.8.0

### Input
nbamfile=$1
tbamfile=$2
prefix=$3
patientID=$4
ILIST=${ILIST:?}
XMX=$6
DEST=$7 ## Destination directory
assert_directory_exists "${DEST}"
XMX=${XMX:-Xmx160G} ## Default 160gb
PADDING=${PADDING:-0} ## Padding the intervals

echo "- clean intermediate files=${CLEAN:?}"

echo "Input:"
echo "- nbamfile=${nbamfile:?}"
echo "- tbamfile=${tbamfile:?}"
echo "- prefix=${prefix:?}"
echo "- patientID=${patientID:?}"
echo "- ILIST=${ILIST:?}"
echo "- PADDING=${PADDING:?}"
echo "- XMX=${XMX:?}"

## Assert existance of input files
assert_file_exists "${nbamfile}"
assert_file_exists "${tbamfile}"
assert_file_exists "${ILIST}"

assert_file_executable "${GATK4:?}"
#assert_file_executable "${LG3_HOME}"/${GATK4}4-funcotator-vcf2tsv

echo "Software:"
python --version
java -version
echo "GATK4 = ${GATK4}"

assert_python ""

### References
assert_file_exists "${REF:?}"
echo "- reference = ${REF}"
assert_file_exists "${GNOMAD}"
echo "- GNOMAD = ${GNOMAD}"

${GATK4} GetSampleName --verbosity ${VERBOSITY} -I "${nbamfile}" -O normal_name.txt
normalname=$(cat normal_name.txt)
echo "Got Normal name: ${normalname}"

${GATK4} GetSampleName --verbosity ${VERBOSITY} -I "${tbamfile}" -O tumor_name.txt
tumorname=$(cat tumor_name.txt)
rm -f tumor_name.txt normal_name.txt
echo "Got Tumor name: ${tumorname}"

echo "-------------------------------------------------"
echo -n "[Mutect2] Somatic Mutation Detection "
date
echo "-------------------------------------------------"
echo "[Mutect2] Patient ID: $patientID"
echo "[Mutect2] Normal bam file: $nbamfile"
echo "[Mutect2] Tumor bam file: $tbamfile"
echo "[Mutect2] Normal Sample: $normalname"
echo "[Mutect2] Tumor Sample: $tumorname"
echo "[Mutect2] Prefix: $prefix"
echo "-------------------------------------------------"
echo "[Mutect2] Java Memory Xmx value: $XMX"
echo -n "[Mutect2] Working directory: "
pwd
echo "-------------------------------------------------"


INT=${ILIST}
echo "Intervals : ${INT}"

echo -e "\\n\\n****** Somatic variations using MuTect2! ******"

{ time ${GATK4} --java-options "-Xms64G -Xmx64G" Mutect2 \
      --verbosity ${VERBOSITY} \
		--tmp-dir "${TMP_DIR}" \
      -R "${REF}" \
      -L "${INT}" \
		-ip "${PADDING}" \
      -I "${tbamfile}" \
      -I "${nbamfile}" \
      -normal "${normalname}" \
      -germline-resource "${GNOMAD}" \
      --f1r2-tar-gz f1r2.tar.gz \
      -O "${prefix}".unfiltered.vcf; } 2>&1 || error "MuTect2 FAILED"
assert_file_exists "${prefix}".unfiltered.vcf
assert_file_exists f1r2.tar.gz
echo "Total calls before filter : "
grep -c "^##" "${prefix}".unfiltered.vcf
echo "****** MuTect2 Completed! ******"


echo -e "\\n\\n****** Pass raw data to LearnReadOrientationModel ******"
{ time ${GATK4} LearnReadOrientationModel \
      --verbosity ${VERBOSITY} \
		--tmp-dir "${TMP_DIR}" \
      -I f1r2.tar.gz \
      -O read-orientation-model.tar.gz; } 2>&1 || error "LearnReadOrientationModel FAILED"
assert_file_exists read-orientation-model.tar.gz
echo "****** LearnReadOrientationModel Completed! ******"

echo -e "\\n\\n****** Normal GetPileupSummaries on known variant sites.******"
{ time ${GATK4} GetPileupSummaries \
      --verbosity ${VERBOSITY} \
		--tmp-dir "${TMP_DIR}" \
      -I "${nbamfile}" \
      -V "${COMMON}" \
      -L "${COMMON}" \
		-ip "${PADDING}" \
      --interval-set-rule INTERSECTION \
      -O "${prefix}".normal_pileup.table; } 2>&1 || error "Normal GetPileupSummaries FAILED"
assert_file_exists "${prefix}".normal_pileup.table
echo "****** Normal GetPileupSummaries Completed! ******"

echo -e "\\n\\n****** Tumor GetPileupSummaries on known variant sites.******"
{ time ${GATK4} GetPileupSummaries \
      --verbosity ${VERBOSITY} \
		--tmp-dir "${TMP_DIR}" \
      -I "${tbamfile}" \
      -V "${COMMON}" \
      -L "${COMMON}" \
		-ip "${PADDING}" \
      --interval-set-rule INTERSECTION \
      -O "${prefix}".tumor_pileup.table; } 2>&1 || error "Tumor GetPileupSummaries FAILED"
assert_file_exists "${prefix}".tumor_pileup.table
echo "****** Tumor GetPileupSummaries Completed! ******"

echo -e "\\n\\n****** Estimate contamination with CalculateContamination.******"
{ time ${GATK4} CalculateContamination \
      --verbosity ${VERBOSITY} \
		--tmp-dir "${TMP_DIR}" \
      -I "${prefix}".tumor_pileup.table \
      --matched-normal "${prefix}".normal_pileup.table \
      --tumor-segmentation "${prefix}".segments.table \
      -O "${prefix}".contamination.table; } 2>&1 || error "CalculateContamination FAILED"
assert_file_exists "${prefix}".segments.table
assert_file_exists "${prefix}".contamination.table
echo "Contamination discovered: "
cat "${prefix}".contamination.table
echo "****** CalculateContamination Completed! ******"

echo -e "\\n\\n****** Pass learned read orientation model to FilterMutectCallswith ******"
{ time ${GATK4} FilterMutectCalls \
      --verbosity ${VERBOSITY} \
		--tmp-dir "${TMP_DIR}" \
      -R "${REF}" \
      -L "${INT}" \
		-ip "${PADDING}" \
      -V "${prefix}".unfiltered.vcf \
      --tumor-segmentation "${prefix}".segments.table \
      --contamination-table "${prefix}".contamination.table \
      --ob-priors read-orientation-model.tar.gz \
      -O "${prefix}".filtered.vcf;  } 2>&1 || error "FilterMutectCalls FAILED"
assert_file_exists "${prefix}".filtered.vcf
echo -ne "Total filtered/PASSed calls:  "
grep -c PASS "${prefix}".filtered.vcf
echo "****** FilterMutectCalls Completed! ******"

### Somatic data source for Funcotator
assert_directory_exists "${FUNCOTATOR:?}"
echo "- FUNCO_PATH=${FUNCOTATOR}"

echo -e "\\n\\n****** Run Funcotator ******"
{ time ${GATK4} Funcotator \
      --verbosity ${VERBOSITY} \
		--tmp-dir "${TMP_DIR}" \
      --variant "${prefix}".filtered.vcf \
      --reference "${REF}" \
      --ref-version hg38 \
      --data-sources-path "${FUNCOTATOR}" \
      --output "${prefix}".variants.funcotated.vcf \
      --output-file-format VCF; } 2>&1 || error "Funcotator FAILED"
assert_file_exists "${prefix}".variants.funcotated.vcf
echo -n "Total PASSed and funcotated calls : "
grep -c PASS "${prefix}".variants.funcotated.vcf 
echo "****** Funcotator Completed! ******"

echo -e "\\n\\n****** Extracting Funcotator annotations in .tsv format ******"
OUT="${prefix}".variants.funcotated.tsv
#"${LG3_HOME}"/${GATK4}4-funcotator-vcf2tsv "${OUT}"

{ time ${GATK4} VariantsToTable \
      --verbosity ${VERBOSITY} \
		--tmp-dir "${TMP_DIR}" \
      -R "${REF}" \
      -V "${prefix}".variants.funcotated.vcf \
      -F FUNCOTATION \
      -O "${OUT}"; } 2>&1 || { echo "FAILED"; exit 1; }
echo -n "Output : "
wc -l "${OUT}"
sed -i 's/|/\t/g;s/\],\[/\n/g;s/\[//;s/\]//' "${OUT}"

NG=$(tail -n +2 "${OUT}" | cut -f1 | sort | uniq -c | wc -l)
echo "List of mutated genes : ${NG}"
tail -n +2 "${OUT}" | cut -f1 | sort | uniq -c

echo "****** VariantsToTable Completed! ******"

if ${CLEAN}; then
   echo "****** Cleaning ... *******"
   rm -f f1r2.tar.gz
   rm -f read-orientation-model.tar.gz
   rm -f "${prefix}".normal_pileup.table
   rm -f "${prefix}".tumor_pileup.table
   rm -f "${prefix}".unfiltered.vcf
   rm -f "${prefix}".unfiltered.vcf.idx
   rm -f "${prefix}".unfiltered.vcf.stats
   rm -f "${prefix}".segments.table
   echo "****** Cleaning Complete! *******"
fi

echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] END: $PROGRAM"
