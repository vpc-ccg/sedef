# SEDEF: SEgmental Duplication Evaluation Framework

# UNDER CONSTRUCTION

## Compile

Simple! Use

```bash
git clone https://github.com/vpc-ccg/sedef
cd sedef
make -j release
```

By default, SEDEF uses Intel C++ compiler. If you are using g++, build with:

```bash
make -j release CXX=g++
```

## Run

Suppose that our genome is in `hg19.fa` file.

First index the file:

```
samtools faidx hg19.fa
```

Then run the `sedef-search` in parallel to get the initial hits:

```bash
# VPC
#echo "qsub -cwd -V -b y -N \"S_${i}_${j}_${m}\" -l h_vmem=10G -l h_rt=24:00:00 -l h_stack=8M " \
#	"python2.7 mesa sedef/sedef search single fasta/hg19.fa chr$i chr$j $m"
for i in `seq 1 22` X Y; do 
for j in `seq 1 22` X Y; do  
	SI=`awk '$1=="chr'$i'" {print $2}' hg19.fa.fai`; 
	SJ=`awk '$1=="chr'$j'" {print $2}' hg19.fa.fai`; 
	if [ "$SI" -le "$SJ" ] ; 
	then 
		for m in y n ; do
		echo "~/mesa ./sedef search single hg19.fa chr$i chr$j $m >out/${i}_${j}_${m}.bed 2>out/log/${i}_${j}_${m}.log"
		done; 
	fi
done
done | time parallel --will-cite -j 80 --eta
>> 10m 11s

grep Total out/log/*.log | wc -l
>> 600s
grep Wall out/log/*.log | tr -d '(' | awk '{s+=$4}END{print s}'
>> 43414.7 (12.06 h)
grep Memory out/log/*.log | awk '{if($3>m)m=$3}END{print m}'
>> 4889.0
```

Then use `sedef-align` to bucket the files for the optimal parallel alignment, and
afterwards run the whole alignment:

```bash
~/mesa ./sedef align bucket out out/bins 1000
>> 36s (5.2G)

for j in out/bins/bucket_???? ; do
	k=$(basename $j);
	echo "~/mesa ./sedef align generate hg19.fa $j 11 >${j}.bed 2>out/log/bins/${k}.log"
done | time parallel --will-cite -j 80 --eta
>> 5:06
503
grep Finished out/log/bins/*.log | wc -l
>> 1000
grep Wall out/log/bins/*.log | tr -d '(' | awk '{s+=$4}END{print s}'
>> 16919.4 (4.70)
17086
grep Memory out/log/bins/*.log | awk '{if($3>m)m=$3}END{print m}'
>> 4993.0
```

Finally, run `sedef-stats` to produce the final output:

```bash
cat out/*.bed > out.bed
cat out/bins/bucket_???? > out.init.bed
cat out/bins/*.bed > out.final.bed
wc -l out.*bed
>>  1656305 out/out.bed
>>  1558896 out/out.init.bed
>>   231472 out/out.final.bed
>>	  67467 out.hg19.bed
```

Then analyse:

```bash
rsync -Pva meganode:/local-scratch/ibrahim/sedef/out.*bed  out/
for i in out/out.bed out/out.init.bed out/out.final.bed ; do 
	mesa python scratch/check-overlap.py $i > ${i}.log  ; 
done

~/mesa ./sedef stats hg19.fa out.final.bed > out.hg19.bed
>> 0m 58s (7.5G)
```

