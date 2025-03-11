SAMPLES="Z00432.g.vcf Z00433.g.vcf Z00434.g.vcf Z00435.g.vcf Z00436.g.vcf Z00437.g.vcf Z00438.g.vcf"

SAMPLES=$(echo "${SAMPLES[*]}" | sed 's/.g.vcf//g')
echo "${SAMPLES[*]}"

SAMPLES=$(echo "${SAMPLES[*]}" | cut -d. -f1)
echo "${SAMPLES[*]}"


INPUTS_BAM2="/costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/X00217.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/X00218.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/X00219.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/X00220.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/X00221.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/X00222.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/X00223.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/X00224.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/X00225.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/X00226.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/Z00432.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/Z00433.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/Z00434.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/Z00435.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/Z00436.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/Z00437.bwa.mrkDups.sort.recal.bam /costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83/Z00438.bwa.mrkDups.sort.recal.bam"

SAMPLES2=$(basename -a "${INPUTS_BAM2[*]}" | cut -d. -f1)
echo "SAMPLES2 = ${SAMPLES2[*]}"

SAMPLES2=$(basename -a "${INPUTS_BAM2}" | cut -d. -f1)
echo "SAMPLES2 = ${SAMPLES2}"
