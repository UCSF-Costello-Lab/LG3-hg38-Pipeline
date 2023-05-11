#!/bin/bash

# shellcheck source=scripts/utils.sh
source "${LG3_HOME:?}/scripts/utils.sh"
source_lg3_conf
XMX=${XMX:-Xmx128G} 

PROGRAM=${BASH_SOURCE[0]}
echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] BEGIN: $PROGRAM"
echo "Call: ${BASH_SOURCE[*]}"
echo "Script: $PROGRAM"
echo "Arguments: $*"
CLEAN=true
VERBOSITY=ERROR ## WARNING INFO DEBUG

### Configuration
LG3_HOME=${LG3_HOME:?}
LG3_OUTPUT_ROOT=${LG3_OUTPUT_ROOT:-output}
LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT:?}
LG3_DEBUG=${LG3_DEBUG:-true}
ncores=${SLURM_NTASKS:-1}
assert_file_exists "${ILIST:?}"

### Debug
if $LG3_DEBUG ; then
  echo "Settings:"
  echo "- LG3_HOME=$LG3_HOME"
  echo "- LG3_OUTPUT_ROOT=$LG3_OUTPUT_ROOT"
  echo "- LG3_SCRATCH_ROOT=$LG3_SCRATCH_ROOT"
  echo "- PWD=$PWD"
  echo "- USER=$USER"
  echo "- hostname=$(hostname)"
  echo "- ncores=$ncores"
  echo "- XMX=${XMX}"
  echo "- CLEAN=${CLEAN}"
  #echo "- ILIST=${ILIST:?}"
  #echo "- PADDING=${PADDING:?}"
fi

## Positional arguments
fastq1=$1
fastq2=$2
S=$3

OUT=${S}_bwamem.sam

echo "Input:"
echo "- fastq1=${fastq1:?}"
echo "- fastq2=${fastq2:?}"
echo "- SAMPLE=${S:?}"

#module load gatk/4.3.0.0
echo "GATK4 = ${GATK4}"
module load openjdk/1.8.0

assert_file_exists "${REF}"
assert_file_exists "${KNOWN_SNPS}"
assert_file_exists "${KNOWN_INDELS}"
assert_file_exists "${KNOWN_INDELS2}"

assert_file_exists ${fastq1}
assert_file_exists ${fastq2}


echo -e "\\n\\n****** BWA-MEM Alignment ******"
## -M ##  flag shorter split hits as secondary (for Picard)
## -p ## smart pairing (ignoring in2.fq)
## -Y ## Use soft clipping for supplementary alignments
## -K 100000000 ## process INT input bases in each batch regardless of nThreads (for reproducibility)
{ time bwa mem -v3 -M -t ${ncores} -K 100000000 \
   "${REF}" "${fastq1}" "${fastq2}" > "${OUT}"; } 2>&1 || error "BWA-MEM FAILED"
assert_file_exists ${OUT}

#samtools flagstat ${OUT} > ${OUT}.flagstat
#cat ${OUT}.flagstat
echo "****** BWA-MEM Completed! ******"


echo -e "\\n\\n****** RevertSam : Creating unmapped uBAM from ${S}_bwamem.sam ******"
{ time ${GATK4} RevertSam \
   -I ${S}_bwamem.sam \
   -O ${S}_u.bam \
   -ATTRIBUTE_TO_CLEAR XS -ATTRIBUTE_TO_CLEAR XA; } 2>&1 || error "RevertSam FAILED"
assert_file_exists ${S}_u.bam
echo "****** RevertSam Completed! ******"

echo -e "\\n\\n****** AddOrReplaceReadGroups : Add read group to ${S}_u.bam ******"
{ time ${GATK4} AddOrReplaceReadGroups \
   -I ${S}_u.bam \
   -O ${S}_rg.bam \
   --RGID ${S} \
   --RGSM ${S} \
   --RGLB ${S} \
   --RGPU Exome \
   --RGPL Illumina; } 2>&1 || error "AddOrReplaceReadGroups FAILED"
assert_file_exists ${S}_rg.bam
echo "****** AddOrReplaceReadGroups Completed! ******"


echo -e "\\n\\n****** MergeBamAlignment : Merge uBAM with aligned BAM ******"
{ time ${GATK4} MergeBamAlignment \
   -ALIGNED ${S}_bwamem.sam \
   -UNMAPPED ${S}_rg.bam \
   -O ${S}_m.bam \
   -R "${REF}" \
   --SORT_ORDER unsorted \
   --CLIP_ADAPTERS false \
   --ADD_MATE_CIGAR true \
   --MAX_INSERTIONS_OR_DELETIONS -1 \
   --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
   --UNMAP_CONTAMINANT_READS false \
   --MIN_UNCLIPPED_BASES 32 \
   --ATTRIBUTES_TO_RETAIN XS -ATTRIBUTES_TO_RETAIN XA; } 2>&1 || error "MergeBamAlignment FAILED"
   assert_file_exists ${S}_m.bam
echo "****** MergeBamAlignment Completed! ******"


echo -e "\\n\\n****** MarkDuplicates : Flag duplicate reads ******"
{ time ${GATK4} MarkDuplicates \
   -I ${S}_m.bam \
   -O ${S}_md.bam \
   --METRICS_FILE ${S}_dup_metrics.txt \
   --REMOVE_DUPLICATES false \
   --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
   --ASSUME_SORT_ORDER queryname; } 2>&1 || error "MarkDuplicates FAILED"
assert_file_exists ${S}_md.bam

samtools flagstat ${S}_md.bam > ${S}_md.bam.flagstat
echo "Flagstat after MarkDuplicates"
cat ${S}_md.bam.flagstat

echo "****** MarkDuplicates Completed! ******"

## Setup pipe
set -o pipefail
echo -e "\\n\\n****** SortSam|SetNmMdAndUqTags : Sort, fix tags and index clean BAM ******"
{ time ${GATK4} SortSam \
   --VERBOSITY ${VERBOSITY} \
   -I ${S}_md.bam \
   -O /dev/stdout \
   --SORT_ORDER coordinate | \
   ${GATK4} SetNmMdAndUqTags \
   --VERBOSITY ${VERBOSITY} \
   -I /dev/stdin \
   -O ${S}_snaut.bam \
   --CREATE_INDEX true \
   -R "${REF}"; } 2>&1 || error "SortSam|SetNmMdAndUqTags FAILED"
assert_file_exists ${S}_snaut.bam
echo "****** SortSam|SetNmMdAndUqTags Completed! ******"


echo -e "\\n\\n****** BaseRecalibrator: first pass ******"
{ time ${GATK4} BaseRecalibrator \
   --verbosity ${VERBOSITY} \
   -I ${S}_snaut.bam \
   -R "${REF}" \
   --known-sites "${KNOWN_SNPS}" \
   --known-sites "${KNOWN_INDELS}" \
   --known-sites "${KNOWN_INDELS2}" \
   -O ${S}_recal_data.table; } 2>&1 || error "BaseRecalibrator FAILED"
assert_file_exists ${S}_recal_data.table
echo "****** BaseRecalibrator Completed! ******"

echo -e "\\n\\n****** ApplyBQSR : apply recalibration model ******"
{ time ${GATK4} ApplyBQSR \
   --verbosity ${VERBOSITY} \
   -I ${S}_snaut.bam \
   --bqsr-recal-file ${S}_recal_data.table \
   -O ${S}.${RECAL_BAM_EXT}.bam; } 2>&1 || error "ApplyBQSR FAILED"
assert_file_exists ${S}.bwa.mrkDups.sort.recal.bam
echo "****** ApplyBQSR Completed! ******"

echo -e "\\n\\n****** ValidateSamFile : Validate final BAM ******"
{ time ${GATK4} ValidateSamFile \
   --VERBOSITY ${VERBOSITY} \
   -R ${REF} \
   -I ${S}.${RECAL_BAM_EXT}.bam; } 2>&1 || error "ValidateSamFile FAILED"
echo "****** ValidateSamFile Completed! ******"

echo -e "\\n\\n****** CollectMultipleMetrics : ******"
{ time ${GATK4} CollectMultipleMetrics \
   --VERBOSITY ${VERBOSITY} \
   -R ${REF} \
	-O ${S}_multi \
   -I ${S}.${RECAL_BAM_EXT}.bam; } 2>&1 || error "CollectMultipleMetrics FAILED"
echo "****** CollectMultipleMetrics Completed! ******"

## Included in CollectMultipleMetrics!
#echo -e "\\n\\n****** CollectInsertSizeMetrics : ******"
#{ time ${GATK4} CollectInsertSizeMetrics \
   #--VERBOSITY ${VERBOSITY} \
   #-R ${REF} \
   #-O ${S}_insert_size_metrics.txt \
	#-H ${S}_insert_size_metrics.pdf \
   #-I ${S}.bwa.mrkDups.sort.recal.bam; } 2>&1 || error "CollectInsertSizeMetrics FAILED"
#echo "****** CollectInsertSizeMetrics Completed! ******"

## Require interval list with header!
echo -e "\\n\\n****** CollectHsMetrics : ******"
{ time ${GATK4} CollectHsMetrics \
   --VERBOSITY ${VERBOSITY} \
   -R ${REF} \
   -O ${S}_Hs_metrics.txt \
	--BAIT_INTERVALS ${ILIST2} \
	--TARGET_INTERVALS ${ILIST2} \
   -I ${S}.${RECAL_BAM_EXT}.bam; } 2>&1 || error "CollectHsMetrics FAILED"
echo "****** CollectHsMetrics Completed! ******"

samtools flagstat ${S}.${RECAL_BAM_EXT}.bam > ${S}_final.bam.flagstat
echo "Flagstat after MarkDuplicates"
cat ${S}_final.bam.flagstat

if $CLEAN ; then
	echo "Deleting tmp files ..."
	rm -rf ${S}_bwamem.sam
	rm -rf ${S}_u.bam
	rm -rf ${S}_rg.bam
	rm -rf ${S}_m.bam
	rm -rf ${S}_md.bam
	rm -rf ${S}_snaut.bam
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] END: $PROGRAM"