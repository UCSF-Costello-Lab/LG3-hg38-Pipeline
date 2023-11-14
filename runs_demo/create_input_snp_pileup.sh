#!/bin/bash

echo "usage: %s patient_ID_conversion.txt bampath"

python /c4/home/jocostello/repos/LG3-hg38-Pipeline/scripts/create_input_snp_pileup.py \
		patient_ID_conversions_P83_Z.txt \
		/costellolab/data4/LG3_hg38_output/LG4/bwa-mem/Patient83

