# SEDEF: SEgmental Duplication Evaluation Framework

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
done | parallel --will-cite -j 80 --eta
>> 10m 32s

grep Total out/log/*.log | wc -l
>> 600s
grep Wall out/log/*.log | tr -d '(' | awk '{s+=$4}END{print s}'
>> 44880s (12.47 h)
```

Then use `sedef-align` to bucket the files for the optimal parallel alignment, and
afterwards run the whole alignment:

```bash
~/mesa ./sedef align bucket out out/bins 1000
>> 47s

for j in out/bins/bucket_???? ; do
	k=$(basename $j);
	echo "~/mesa ./sedef align generate hg19.fa $j 11 >${j}.bed 2>out/log/bins/${k}.log"
done | time parallel --will-cite -j 80 --eta
>> 27m 33s

grep Finished out/log/bins/*.log | wc -l
>> 1000
grep Wall out/log/bins/*.log | tr -d '(' | awk '{s+=$4}END{print s}'
>> 43704 (12.14 h)  
```

Finally, run `sedef-stats` to produce the final output:

```bash
cat out/bins/bucket_???? > out.init.bed
cat out/bins/*.bed > out.final.bed
wc -l out.*bed
>>   975511 out.final.bed
>>  1656305 out.init.bed
```

### Mouse

```bash
for i in `seq 1 19` X Y; do 
	for j in `seq 1 19` X Y; do  
		SI=`awk '$1=="chr'$i'" {print $2}' mm8.fa.fai`; 
		SJ=`awk '$1=="chr'$j'" {print $2}' mm8.fa.fai`; 
		if [ "$SI" -le "$SJ" ] ; then 
			for m in y n ; do
			echo "~/mesa ./sedef search single mm8.fa chr$i chr$j $m >mouse_out/${i}_${j}_${m}.bed 2>mouse_out/log/${i}_${j}_${m}.log"
			done; 
		fi
	done
done | time parallel --will-cite -j 80 --eta
>> 12m 17s
grep Total mouse_out/log/*.log | wc -l
>> 462
grep Wall mouse_out/log/*.log | tr -d '(' | awk '{s+=$4}END{print s}'
>> 46962.9s (13.05 h)
~/mesa ./sedef align bucket mouse_out mouse_out/bins 1000
>> 
for j in mouse_out/bins/bucket_???? ; do
	k=$(basename $j);
	echo "~/mesa ./sedef align generate mm8.fa $j 11 >${j}.bed 2>mouse_out/log/bins/${k}.log"
done | time parallel --will-cite -j 80 --eta
>> 
grep Finished mouse_out/log/bins/*.log | wc -l
>> 1000
grep Wall mouse_out/log/bins/*.log | tr -d '(' | awk '{s+=$4}END{print s}'
>> 236095 (65.58 h)cv
cat mouse_out/bins/bucket_???? > mouse_out.init.bed
cat mouse_out/bins/*.bed > mouse_out.final.bed
wc -l mouse_out.*bed
>>   975511 mouse_out.final.bed
>>  1656305 mouse_out.init.bed

# Merge logs and remove progress bars and headers

for S in G S X Y ; do
	for i in ${S}_*.e*; do grep -v '%' $i | tail -n+4 ; done > ${S}.log
	cat ${S}_*.o* > ${S}.bed
	rm -rf ${S}*.[eo]*
done

for i in G X Y S ; do 
	mkdir -p bins/${i}
	sedef/sedef align bucket ${i}.bed bins/${i} 2000 
done

for i in G X Y S ; do 
	for j in bins/${i}/* ; do
		k=$(basename $j);
		echo "qsub -cwd -V -b y -N \"${i}_${k}\" -l h_vmem=10G -l h_rt=24:00:00 -l h_stack=8M " \
			"python2.7 mesa sedef/sedef align generate fasta/hg19.fa $j"
	done
done

for S in G X Y S ; do 
	for i in ${S}_*.e*; do grep -v '\.\.' $i | tail -n+4 ; done > ${S}.align.log
	cat ${S}_*.o* > ${S}.align.bed
	rm -rf ${S}_*.[eo]*
done

# Count stuff
for S in G S X Y ; do 
	SEA=`cat ${S}.log | grep 'Wall' | awk '{print $4}' | tr -d '[()]' | awk '{s+=$1} END{print s}'`
	ALN=`cat ${S}.align.log | grep 'Wall' | awk '{print $4}' | tr -d '[()]' | awk '{s+=$1} END{print s}'`
	printf "%s: %5.1f %5.1f %5.1f\n" $S $((SEA/3600)) $((ALN/3600)) $(((SEA+ALN)/3600))
done

#mkdir -p output/search
#zmv '(sedefrun_*).o*' 'output/search/$1.bed'
```


