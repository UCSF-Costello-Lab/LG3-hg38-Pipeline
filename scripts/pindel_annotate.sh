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
CLEAN=false

### Debug
if $LG3_DEBUG ; then
  echo "Settings:"
  echo "- LG3_HOME=$LG3_HOME"
  echo "- LG3_OUTPUT_ROOT=$LG3_OUTPUT_ROOT"
  echo "- LG3_SCRATCH_ROOT=$LG3_SCRATCH_ROOT"
  echo "- PWD=$PWD"
  echo "- USER=$USER"
fi

PROG=$(basename "$0")
OK() {
        echo "OK: line ${BASH_LINENO[0]} in $PROG"
}

datafile=$1
proj=$2
echo "Input:"
echo "- datafile=${datafile:?}"
echo "- proj=${proj:?}"

assert_file_exists "${datafile}"

BIN=${LG3_HOME}/scripts
assert_directory_exists "${BIN}"
assert_directory_exists "${ANNOVAR_HOME}"

#KINASEDATA="${LG3_HOME}/resources/all_human_kinases.txt"
#COSMICDATA="${LG3_HOME}/resources/CosmicMutantExport_v58_150312.tsv"
#CANCERDATA="${LG3_HOME}/resources/SangerCancerGeneCensus_2012-03-15.txt"
#CONVERT="${LG3_HOME}/resources/RefSeq.Entrez.hg38.2023-05-23.txt
ANNDB=${ANNOVAR_HOME}/hg38db/

echo "References:"
echo "- KINASEDATA=${KINASEDATA:?}"
echo "- COSMICDATA=${COSMICDATA:?}"
echo "- CANCERDATA=${CANCERDATA:?}"
echo "- CONVERT=${CONVERT:?}"
echo "- ANNDB=${ANNDB:?}"

assert_file_exists "${KINASEDATA}"
assert_file_exists "${COSMICDATA}"
assert_file_exists "${CANCERDATA}"
assert_file_exists "${CONVERT}"
assert_directory_exists "${ANNDB}"

echo -n "Started $PROG on "
date

### run AnnoVar
echo "================= [Annotate] run annovar"
"$ANNOVAR_HOME/annotate_variation.pl" -filter -dbtype 1000g2014oct_all -buildver hg38 "${datafile}.filter.intersect" "$ANNDB" || error "annotate_variation.pl failed"
OK 

awk -F '\t' '{for(i=3;i<=NF;i++) {printf $i"\t";} print $1}' "${datafile}.filter.intersect.hg38_ALL.sites.2014_10_dropped" > "${datafile}.tmp11"
assert_file_exists "${datafile}.tmp11"

awk -F '\t' '{for(i=1;i<=NF;i++) {printf $i"\t";} print ""}' "${datafile}.filter.intersect.hg38_ALL.sites.2014_10_filtered" > "${datafile}.tmp12"
assert_file_exists "${datafile}.tmp12"

cat "${datafile}.tmp11" "${datafile}.tmp12" > "${datafile}.tmp1"
assert_file_exists "${datafile}.tmp1"

"$ANNOVAR_HOME/annotate_variation.pl" -filter -dbtype 1000g2015aug_all -buildver hg38 "${datafile}.tmp1" "$ANNDB" || error "annotate_variation.pl failed"
OK 

awk -F '\t' '{for(i=3;i<=NF;i++) {printf $i"\t";} print $1}' "${datafile}.tmp1.hg38_ALL.sites.2015_08_dropped" > "${datafile}.tmp21"
assert_file_exists "${datafile}.tmp21"

awk -F '\t' '{for(i=1;i<=NF;i++) {printf $i"\t";} print ""}' "${datafile}.tmp1.hg38_ALL.sites.2015_08_filtered" > "${datafile}.tmp22"
assert_file_exists "${datafile}.tmp22"

cat "${datafile}.tmp21" "${datafile}.tmp22" > "${datafile}.tmp2"
assert_file_exists "${datafile}.tmp2"

"$ANNOVAR_HOME/annotate_variation.pl" -filter -dbtype avsnp150 -buildver hg38 "${datafile}.tmp2" "$ANNDB" || error "annotate_variation.pl failed"
OK 

awk -F '\t' '{for(i=3;i<=NF;i++) {printf $i"\t";} print $1"_"$2}' "${datafile}.tmp2.hg38_avsnp150_dropped" > "${datafile}.tmp31"
assert_file_exists "${datafile}.tmp31"

awk -F '\t' '{for(i=1;i<=NF;i++) {printf $i"\t";} print ""}' "${datafile}.tmp2.hg38_avsnp150_filtered" > "${datafile}.tmp32"
assert_file_exists "${datafile}.tmp32"

cat "${datafile}.tmp31" "${datafile}.tmp32" > "${datafile}.tmp3"
assert_file_exists "${datafile}.tmp3"

"$ANNOVAR_HOME/annotate_variation.pl" --geneanno --buildver hg38 --outfile "${datafile}.filter.intersect.anno" "${datafile}.tmp3" "$ANNDB" || error "annotate_variation.pl failed"
OK 

### clean up AnnoVar exonic data, put into final mutation table format
echo "================= [Annotate] reformat annovar"
"$BIN/pindel_reformat_annovar.py" "${datafile}.filter.intersect.anno" "${datafile}.filter" || error "pindel_reformat_annovar.py failed"
assert_file_exists "${datafile}.filter"
OK 

### annotate with normal coverage
echo "================= [Annotate] annotate with normal coverage"
"$BIN/pindel_annotate_normal_coverage.py" "${datafile}.filter.intersect.anno.muts" "${proj}" || error "pindel_annotate_normal_coverage.py failed"
assert_file_exists "${datafile}.filter.intersect.anno.muts"
OK 

## annotate with Kinase & Cosmic & Sanger Cancer Gene
echo "================= [Annotate] annotate with cosmic, kinase, sanger cancer gene list"
"$BIN/annotation_COSMIC.py" "${datafile}.filter.intersect.anno.muts.norm.txt" "$COSMICDATA" > "${datafile}.filter.intersect.anno.muts.tmp1" || error "annotation_COSMIC.py failed"
assert_file_exists "${datafile}.filter.intersect.anno.muts.tmp1"
OK 

"$BIN/annotation_KINASE.py" "${datafile}.filter.intersect.anno.muts.tmp1" "$KINASEDATA" > "${datafile}.filter.intersect.anno.muts.tmp2" || error "annotation_KINASE.py failed"
assert_file_exists "${datafile}.filter.intersect.anno.muts.tmp2"
OK 

"$BIN/annotation_CANCER.py"  "${datafile}.filter.intersect.anno.muts.tmp2" "$CANCERDATA" "$CONVERT" >  "${datafile}.filter.intersect.anno.muts.norm.anno.txt" || error "annotation_CANCER.py failed"
assert_file_exists "${datafile}.filter.intersect.anno.muts.norm.anno.txt"
OK 

## remove indels with <14 reads of raw coverage in the normal
echo "================= [Annotate] remove indels with <14 reads in normal"
awk -F'\t' '{if($21>=14) print}' "${datafile}.filter.intersect.anno.muts.norm.anno.txt" > "${datafile}.filtered.anno.txt"
assert_file_exists "${datafile}.filtered.anno.txt"
OK

if ${CLEAN}; then
### clean up intermediate files
	echo "================= [Annotate] delete intermediate files"
	rm -f "${datafile}.filter"
	rm -f "${datafile}.filter".*
	rm -f "${datafile}.tmp"*
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S %Z')] END: $PROGRAM"
