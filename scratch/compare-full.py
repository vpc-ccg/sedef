#!/bin/bash

WGAC="$1"
SEDEF="$2"

cat $1 \
	| awk '{print $1,$2,$3,$7,$8,$9,$4,$5,"+",$6}' 'OFS=\t' \
	| grep -v '_gl' | sed 's/_$/-/g' | sed 1d > "${WGAC}.bed"

echo "Bedtools..."
bedtools pairtopair -a ${SEDEF} -b "${WGAC}.bed" \
	| cut -f16- | sort -u > "${SEDEF}.wgac_matches.bed" 

comm -23 <(sort -u "${WGAC}.bed") "${SEDEF}.wgac_matches.bed" > "${SEDEF}.wgac_misses.bed" 

echo -ne "Total WGAC hits: "
cat "${WGAC}.bed" | awk '{print $1":"$2"\t"$4":"$5}' \
	| awk '{if($1<$2)print $1,$2; else print $2,$1;}' | sort -u | wc -l 

echo -ne "Total misses: "
cat "${SEDEF}.wgac_misses.bed"  \
	| awk '{print $1":"$2"\\\\|"$4":"$5}' \
	| while read i
	  do 
		grep -q "$i" "${SEDEF}.wgac_matches.bed" || echo $i 
	  done \
	| sed 's/\\|/\t/' | awk '{if($1<$2)print $1,$2; else print $2,$1;}' \
	| sort -u | tee "${SEDEF}.wgac_misses.names" | wc -l

