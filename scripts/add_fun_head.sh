#!/bin/bash

if [ $# -lt 1 ]; then
	echo "Usage: $0 *.variants.funcotated.tsv"	
	exit 1
fi

for f in $@
do
	prefix=$(basename $f .tsv)
	echo "Adding headers to ${prefix}"
	VCF=${prefix}.vcf
	if [ ! -f "${VCF}" ]; then
   	echo "ERROR: Can't find ${VCF}"
   	exit 1
	else
		echo "Found ${VCF}"
	fi
	TSV=${prefix}.tsv
	if [ ! -f "${TSV}" ]; then
   	echo "ERROR: Can't find ${TSV}"
   	exit 1
	else
		echo -n "Found total muts: "
		tail -n +2 "${TSV}" | wc -l
	fi
	
	doit() {
		while read -r N ID
		do
			echo -n "$N -- $ID == "
			cut -f${N} ${TSV} | grep -v FUNCOTATION | grep -cve '^\s*$'
		done
	}
	HEAD=$(grep "ID=FUNCOTATION" ${VCF} | cut -c1-130 --complement )
	echo ${HEAD%??} | tr '|' '\n' | nl | doit > ${prefix}.head
	wc -l ${prefix}.head
	
	echo ${HEAD%??} | tr '|' '\t' > ${prefix}.tab
	echo -n "First/Last columns = " 
	cut -f1,165 ${prefix}.tab
	
	tail -n +2 ${prefix}.tsv >> ${prefix}.tab
	mv ${prefix}.tab ${prefix}.tsv
	wc -l ${prefix}.tsv
	
done

