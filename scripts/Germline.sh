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

### Debug
if $LG3_DEBUG ; then
  echo "$PROG Settings:"
  echo "- LG3_HOME=$LG3_HOME"
  echo "- LG3_OUTPUT_ROOT=$LG3_OUTPUT_ROOT"
  echo "- LG3_SCRATCH_ROOT=$LG3_SCRATCH_ROOT"
  echo "- PWD=$PWD"
  echo "- USER=$USER"
  echo "- hostname=$(hostname)"
  echo "- ncores=$ncores"
fi

## Input
nbamfile=$1
PATIENT=$2
ILIST=$3
echo "Input:"
echo "- nbamfile=${nbamfile:?}"
echo "- PATIENT=${PATIENT:?}"
echo "- ILIST=${ILIST:?}"

## Assert existance of input files
assert_file_exists "${nbamfile}"
assert_file_exists "${ILIST}"

bamdir=${nbamfile%/*}

${GATK4} GetSampleName --verbosity ${VERBOSITY} -I "${nbamfile}" -O normal_name.txt
normalname=$(cat normal_name.txt)
echo "Got Normal name: ${normalname}"
rm normal_name.txt

## References
echo "References:"
echo "- REF=${REF:?}"
#echo "- DBSNP=${DBSNP:?}"
assert_file_exists "${REF}"
#assert_file_exists "${DBSNP}"

## Software
PYTHON_VCF_GERMLINE=${LG3_HOME}/scripts/vcf_germline.py
echo "Software:"
#echo "- JAVA=${JAVA:?}"
#echo "- PYTHON=${PYTHON:?}"
echo "- GATK=${GATK4:?}"
#assert_python "$PYTHON"


## Assert existance of software
assert_file_executable "${JAVA}"
assert_file_executable "${PYTHON}"
assert_file_exists "${GATK4}"
assert_file_exists "${PYTHON_VCF_GERMLINE}"

echo "-------------------------------------------------"
echo "[Germline] Germline SNPs and relatedness"
echo "-------------------------------------------------"
echo "[Germline] Patient ID: $PATIENT"
echo "[Germline] Bam file directory: $bamdir"
echo "[Germline] Normal Sample: $normalname"
echo "-------------------------------------------------"

HaplotypeCaller() {
   echo -e "\\n\\n****** HaplotypeCaller : Call SNP and indel variants in ERC mode ******"
	bambase=$(basename "$1")
	ID=${bambase%%.*}
	echo "Processing ${ID} ..."
   { time ${GATK4} HaplotypeCaller \
      --verbosity ${VERBOSITY} \
		-L ${ILIST} \
		-ip ${PADDING} \
      -R "${REF}" \
		-I "$1" \
      -O "${ID}.g.vcf" \
      -ERC GVCF; } 2>&1 || error "HaplotypeCaller FAILED"
   assert_file_exists "${ID}.g.vcf"
	echo -ne "Called ${ID}.g.vcf : "
	grep -cv '^#' ${ID}.g.vcf
   echo "****** HaplotypeCaller Completed! ******"
}
CombineGVCFs() {
	echo -e "\\n\\n****** CombineGVCFs : Combine multi-sample gVCF ******"
	# shellcheck disable=SC2086
	{ time ${GATK4} CombineGVCFs \
   	--verbosity ${VERBOSITY} \
		-L "${ILIST}" \
		-ip ${PADDING} \
   	-R "${REF}" \
   	-O "${PATIENT}".g.vcf \
   	${INPUTS_VCF}; } 2>&1 || error "CombineGVCFs FAILED"
	assert_file_exists "${PATIENT}.g.vcf"
	echo -n "Total combined variants in ${PATIENT}.g.vcf : "
	grep -c '^#' "${PATIENT}.g.vcf"
	echo "****** CombineGVCFs Completed! ******"
}

for bam in "${bamdir}"/*."${RECAL_BAM_EXT}.bam"
do
	echo "BAM = ${bam}"
	HaplotypeCaller "${bam}"
done

## Construct array with one or more '-V <g.vcf>' elements
INPUTS_VCF=()
INPUTS_VCF=$(for i in *.g.vcf; do echo -n " -V $i"; done)
echo "INPUTS_VCF : ${INPUTS_VCF[*]}"

CombineGVCFs 

OUT=${PATIENT}.raw.vcf
echo -e "\\n\\n****** GenotypeGVCFs : Joint genotyping on a single-sample GVCFs ******"
{ time ${GATK4} GenotypeGVCFs \
   --verbosity ${VERBOSITY} \
   -R "${REF}" \
	-L "${ILIST}" \
	-ip ${PADDING} \
   -O "${OUT}" \
   --variant "${PATIENT}.g.vcf"; } 2>&1 || error "GenotypeGVCFs FAILED"
assert_file_exists "${OUT}"
echo -ne "\\nRaw VCF : ${OUT} : "
grep -vc '^#' ${OUT}
echo -e "****** GenotypeGVCFs Completed! ******\\n"


## Construct string with one or more '-I <bam>' elements
INPUTS_BAM=$(for i in "${bamdir}"/*."${RECAL_BAM_EXT}.bam"
do
        assert_file_exists "${i}"
        echo -n "-I $i "
done)
echo "INPUTS_BAM ${INPUTS_BAM}"

echo -e "\\n\\n****** VariantAnnotator : ******"
OUT=${PATIENT}.annotated.vcf
# shellcheck disable=SC2086
	#-L "${PATIENT}.raw.vcf" \
	#-L "${ILIST}" \
	#-ip "${PADDING}" \
{ time ${GATK4} VariantAnnotator \
   --verbosity ${VERBOSITY} \
   -R "${REF}" \
	${INPUTS_BAM} \
	-V "${PATIENT}.raw.vcf" \
	--dbsnp "${KNOWN_SNPS3}" \
	--annotation QualByDepth \
	--annotation RMSMappingQuality \
	--annotation MappingQualityZero \
	--annotation MappingQualityRankSumTest \
	--annotation FisherStrand \
	--annotation HaplotypeFilteringAnnotation \
	--annotation ReadPosRankSumTest \
	--annotation Coverage \
   -O "${OUT}"; } 2>&1 || error "VariantAnnotator FAILED!" 
assert_file_exists "${OUT}"	
echo -ne "\\nAnnotated raw VCF : ${OUT} : "
grep -vc '^#' "${OUT}"
echo -e "****** VariantAnnotator Completed! ******\\n"

echo -e "\\n\\n****** VariantFiltration : ******"
OUT="${PATIENT}.filtered.vcf"
{ time ${GATK4} VariantFiltration \
   --verbosity ${VERBOSITY} \
   -R "${REF}" \
	-V "${PATIENT}.annotated.vcf" \
	-L "${ILIST}" \
	-ip "${PADDING}" \
	--filter-name "QDFilter" \
	--filter-expression "QD < 2.0" \
	--filter-name "MQFilter" \
	--filter-expression "MQ < 40.0" \
	--filter-name "FSFilter" \
	--filter-expression "FS > 60.0" \
	--filter-name "HaplotypeScoreFilter" \
	--filter-expression "HaplotypeScore > 13.0" \
	--filter-name "MQRankSumFilter" \
	--filter-expression "MQRankSum < -12.5" \
	--filter-name "ReadPosFilter" \
	--filter-expression "ReadPosRankSum < -8.0" \
	-O "${OUT}"; } 2>&1 || error "VariantFiltration FAILED!" 
assert_file_exists "${OUT}"
echo -ne "\\nFiltered VCF : ${OUT} : "
grep -vc '^#' "${OUT}"
echo -ne "PASSed variants : "
grep -c PASS "${OUT}" 
echo -e "****** VariantFiltration Completed! ******\\n"


echo -e "\\n\\n****** SelectVariants : SNPs only ******"
OUT="${PATIENT}.filtered.snps.vcf"
{ time ${GATK4} SelectVariants \
   --verbosity ${VERBOSITY} \
   -R "${REF}" \
	-L "${ILIST}" \
	-ip "${PADDING}" \
   -V "${PATIENT}.filtered.vcf" \
	--select-type-to-include SNP \
	-O "${OUT}"; } 2>&1 || error "SelectVariants FAILED!"
assert_file_exists "${OUT}"
echo -ne "\\nFiltered SNPs : ${OUT} : "
grep -vc '^#' "${OUT}"
echo -ne "PASSed SNPs : "
grep -c PASS "${OUT}"
echo -e "****** SelectVariants Completed! ******\\n"


for i in "${bamdir}"/*.bam
do
   bambase=$(basename "${i}")
   tumorname=${bambase%%.*}
   prefix="NOR-${normalname}_vs_${tumorname}"

    echo "[Germline] Checking germline SNPs for sample relatedness: $tumorname vs $normalname"
    $PYTHON "${PYTHON_VCF_GERMLINE}" \
          "${OUT}" \
          "$normalname" \
          "$tumorname" \
       > "${prefix}.germline" || error "Germline analysis failed"
      assert_file_exists "${prefix}.germline"
done
echo "Deleting germline vs. germline comparison..."
rm -f "NOR-${normalname}_vs_${normalname}.germline"


echo "[Germline] Results:"
#grep Tumor -- *.germline
cat *.germline

echo "[Germline] Finished!"
echo "-------------------------------------------------"

if ${CLEAN} ; then
	echo "Deleting temp files..."
	rm -f *.g.vcf
	rm -f *.g.vcf.idx
	rm -f ${PATIENT}.annotated.vcf
	rm -f ${PATIENT}.annotated.vcf.idx
	rm -f ${PATIENT}.raw.vcf
	rm -f ${PATIENT}.raw.vcf.idx
fi


echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] END: $PROGRAM"
