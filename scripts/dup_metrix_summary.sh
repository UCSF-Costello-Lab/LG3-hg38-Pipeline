#!/bin/bash

if [ ! -r "$1" ]; then
	echo "Usage: $0 *_dup_metrics.txt"	
	exit 1
fi

B=$(basename "$1" .txt)
echo ${B}
OUT=${B}_summary

echo "==== Duplication report human readable ====" > ${OUT}
grep -A1 LIBRARY "$1" | trans | head -n 1 >> ${OUT}
grep -A1 LIBRARY "$1" | trans | awk '{ printf("%s\t%'"'"'f\n", $1,$2) }' | tail -n +2 | sed 's/.000000//' >> ${OUT}
head -n 15 ${OUT}
#wc -l ${OUT}
#grep -A1 LIBRARY X00331_dup_metrics.txt | trans > X00331_dup_metrics_summary
