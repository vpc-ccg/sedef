# SEDEF: The Pipeline

Directory structure:

```
+ fasta
	- chr<1..22,X,Y>.<fa|fai>
+ sedef
	* src
	- sedef-jaccard
	- megaalign
+ logs
	- search.log
	- align_time.log
	- align_mem.log
+ output
	+ search
		- search_<chrA>_<chrB>_<y|n>.bed
	+ align
	- compare.py
	- error_rates.py
	- fileIO.py
	+ <sample>
		+ files
			- <tool>_<mod>.<bam|vcf>
			- giab.vcf
		+ phase
			- <tool>_<mod>.<link|unlink|phase>
```

## Get WGAC data

```bash
wget http://humanparalogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab
cat GRCh37GenomicSuperDup.tab \
	| awk '{print $1,$2,$3,$7,$8,$9,$4,$5,"+",$6}' 'OFS=\t' \
	| grep -v '_gl' \
	| sed 's/_$/-/g' | sed 1d > GRCh37GenomicSuperDup.bed
```

## Run SEDEF searchcd se


G S
X Y

```bash

for i in `seq 1 22` X Y; do 
for j in `seq 1 22` X Y; do  
	SI=`wc -c < fasta/chr${i}.fa`; 
	SJ=`wc -c < fasta/chr${j}.fa`; 
	if [ "$SI" -ge "$SJ" ] ; then for m in y n ; do
		echo "qsub -cwd -V -b y -N \"S_${i}_${j}_${m}\" -l h_vmem=10G -l h_rt=24:00:00 -l h_stack=8M " \
			"python2.7 mesa sedef/sedef search single fasta/hg19.fa chr$i chr$j $m"
	done; fi
done; done | parallel

# Merge logs and remove progress bars and headers
for i in sedefrun_*.e*; do grep -v '%' $i | tail -n+4 ; done > logs/search.log

mkdir -p output/search
zmv '(sedefrun_*).o*' 'output/search/$1.bed'
```


### Bucket the calls

```bash
mkdir -p output/search/buckets
rm -rf output/search/buckets/*; 
for i in output/search/*.bed ; do 
	echo $i; 
	awk '{p = int(($3 - $2) / 1000); print $0 >> "output/search/buckets/b_"p;}' $i; 
done

WC=`wc -l output/search/buckets` # then put it to python dict with sublime
./mesa python2.7 eval.py
rm -rf output/search/buckets
```
Total records: 
>> 558,330,096
>> 289,062,876

Total time: 30,231,139
Wall time:     00:13:27 (807.1 s) for _eval.py_ (ram: 8011.9)

```bash
for i in `seq 0 4000`; do 
	echo -ne "$i\t"; 
	for j in align_$i.e* ; do 
		tail -n20 $j | grep '[0-9]s' | tail -n1 | cut -f2 | tr '\n' '\t' | tr -d 's'  
	done
	echo  
done | awk '{w = 0; for (i = 2; i <= NF; i++) w += $i; e="run_"$1; $1 = w; print e, $0; q += w} END {print "Total", q}' OFS='\t'
```

## Check results

```bash
for SEDEF in sedef_500 sedef_500n ; do
	echo ${SEDEF} ;
	bedtools pairtopair -a ${SEDEF}.bed -b GRCh37GenomicSuperDup.bed > matches.${SEDEF}.bed ;
	cat matches.${SEDEF}.bed  | cut -f16- | sort -u > wgac.matches.${SEDEF}.bed ;
	comm -23 <(sort -u GRCh37GenomicSuperDup.bed) wgac.matches.${SEDEF}.bed > wgac.misses.${SEDEF}.bed ;
	cat GRCh37GenomicSuperDup.bed | awk '{print $1":"$2"\t"$4":"$5}'    | awk '{if($1<$2)print $1,$2; else print $2,$1;}' | sort -u | wc -l ;
	cat wgac.misses.${SEDEF}.bed  | awk '{print $1":"$2"\\\\|"$4":"$5}' | while read i; do grep -q "$i" wgac.matches.${SEDEF}.bed || echo $i ; done | sed 's/\\|/\t/' | awk '{if($1<$2)print $1,$2; else print $2,$1;}' | sort -u | tee wgac.misses.${SEDEF}.id | wc -l ;
	echo ;
done
```
24474
65

```bash
sort sedef.bed  -k1,1V -k4,4V -k10,10 -k2,2n -k5,5n -k3,3n -k6,6n > sedef_sorted.bed
```

0.10


for i in ../align/*.bed ; do
bash -c "





for i in `seq 1 22` X Y ; do echo -ne "chr$i\t500000000" >> genome.bed ; done


for i in `seq 4000 -1 0` ; do
	qsub -cwd -V -b y -N align_$i -l h_vmem=6G -l h_rt=24:00:00 -l h_stack=8M \
		./mesa sedef/megaalign fasta/hg19.fa output/search/balanced/bal_$i
done

 python2.7 mesa sedef/sedef-jaccard $i $j $m

pv/pv ok | bedtools pairtopair -a - -b ~/GRCh37GenomicSuperDup.bed > ok_matches

for i in `seq 0 4000`; do echo -ne "$i\t"; for j in align_$i.e* ; do tail -n20 $j | grep '[0-9]s' | tail -n1 | cut -f2 | tr '\n' '\t' | tr -d 's' ; done; echo  ; done | awk '{w=0; for (i=2;i<=NF;i++) w+=$i; e="run_"$1; $1=w; print e,$0; q+=w;} END{print "Total", q}'  OFS='\t' | tee ../logs/align-time.log
cat search.log | grep Wall | perl -pe 's/.+\(([.0-9]+) s.+/\1/g' | awk '{p+=$1} END{print p}'
tail align_*.e* | grep 'Memory' | awk '{if($3>m)m=$3;print} END{print "Max: ", m}' | tee ../logs/align-mem.log



cat missed_sedef.bed  | awk '{print $1":"$2"\\\\|"$4":"$5}' | while read i; do if ! grep -q "$i" shared_wgac.bed; then echo $i; fi ; done   | tee real_misses.id | wc -l




# Run alignment

for i in *.bed ; do job -n "bedtools_$i" -m 2G bash -c "cat $i | sed  's/\t0\t/\t0\t0\t/' | bedtools pairtopair -a $i -b ~/sedef/GRCh37GenomicSuperDup.bed" ; done

cat >eval.sh <<<EOF
>&2 echo "Begin: `date`"
python2.7 eval.py $1 $2 \
	| awk '{print $1,$2,$3,gensub(/\t/, ";", "g", $0),$8; print $4,$5,$6,"B",$9}' 'OFS=\t' \
	| head \
	| bedtools getfasta -name -fi fasta/hg19.fa -bed - -fo - -s \
	| sedef/megaalign
>&2 echo "Done: `date`"
EOF

for i in `seq 0 2000`; do 
	qsub -cwd -V -b y -N align_$i -l h_vmem=4G -l h_rt=24:00:00 -l h_stack=8M bash eval.sh 2000 $i ;
done


# Merge results


