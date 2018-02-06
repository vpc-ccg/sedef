# SEDEF: SEgmental Duplication Evaluation Framework

SEDEF is a tool to find all segmental duplications in the genome.

## Paper Results

[Please find here](results/out.hg19.bed) the calls for hg19.

## How to compile

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

SEDEF requires Boost libraries in order to compile. In case you have non-standard Boost installation, you can still compile as follows:

```bash
CPATH={path_to_boost} make -j release
```

> **Warning:** Only Linux is currently supported; macOS support is planned later on.

## How to run

Suppose that our genome is in `hg19.fa` file (we used UCSC hg19 with "normal" 24 chromosomes without patches (unGl) or random strains (chrXX_random).

### Automatic transmission

Just go to `sedef` directory and run
```bash
./sedef.sh -o <output> -j <jobs> <genome> 
```

For example, to run hg19.fa on 80 cores type:
```bash
./sedef.sh -o sedef_hg19 -j 80 hg19.fa 
```

You can add `-f` if `sedef_hg19` already exists (it will overwrite its content though). The final results will be
located in `sedef_hg19/final.bed`.

Please note that `sedef.sh` requires SAMtools and GNU Parallel.

### Manual transmission

First make sure to index the file:

```bash
samtools faidx hg19.fa
```

Then run the `sedef-search` in parallel (in this example, we use GNU parallel) to get the initial SD seeds:
```bash
mkdir -p out # For the output
mkdir -p out/log # For the logs

for i in `seq 1 22` X Y; do 
for j in `seq 1 22` X Y; do  
	SI=`awk '$1=="chr'$i'" {print $2}' hg19.fa.fai`; 
	SJ=`awk '$1=="chr'$j'" {print $2}' hg19.fa.fai`; 
	if [ "$SI" -le "$SJ" ] ; 
	then 
		for m in y n ; do
		echo "sedef search single hg19.fa chr$i chr$j $m >out/${i}_${j}_${m}.bed 2>out/log/${i}_${j}_${m}.log"
		done; 
	fi
done
done | time parallel --will-cite -j 80 --eta
# Running time in our case: 10m 11s

# Now make sure that all runs completed successfully 
grep Total out/log/*.log | wc -l
# You should see here 600 (or n(n+1) if you have n chromosomes in your file)

# Get the single-core running time
grep Wall out/log/*.log | tr -d '(' | awk '{s+=$4}END{print s}'
# In our case we get 43414.7 (12.06 h)

# Get the maximum meory usage as well
grep Memory out/log/*.log | awk '{if($3>m)m=$3}END{print m}'
# We got 4889.0
```

Then use `sedef-align` to bucket the files for the optimal parallel alignment, and
afterwards run the whole alignment:
```bash
# First bucket the reads into 1000 bins
mkdir -p out/bins
mkdir -p out/log/bins
time sedef align bucket out out/bins 1000
# Takes 36s (uses 5.2G of memory)

# Now run the alignment
for j in out/bins/bucket_???? ; do
	k=$(basename $j);
	echo "sedef align generate hg19.fa $j 11 >${j}.bed 2>out/log/bins/${k}.log"
done | time parallel --will-cite -j 80 --eta
# It took 5m 06s

# Make sure that all runs finished nicely
grep Finished out/log/bins/*.log | wc -l
# Should be number of bins (in our case, 1000)

# Get again the total running time
grep Wall out/log/bins/*.log | tr -d '(' | awk '{s+=$4}END{print s}'
# We have 17086 (4.74h)

# And the memory
grep Memory out/log/bins/*.log | awk '{if($3>m)m=$3}END{print m}'
# We got 5017.5
```

Finally, run `sedef-stats` to produce the final output:
```bash
# Concatenate the files 
cat out/*.bed > out.bed # seed SDs
cat out/bins/bucket_???? > out.init.bed # potential SD regions
cat out/bins/*.bed > out.final.bed # final chains

# Count the number of SDs in each stage
wc -l out.*bed
#  1656305 out/out.bed
#  1558896 out/out.init.bed
#   231472 out/out.final.bed

# Now get the final calls
sedef stats generate hg19.fa out.final.bed > out.hg19.bed
# Took 0m 58s (7.5G)
```

Final calls will be in `out.hg19.bed`.

## Notes ... 

```# VPC
#echo "qsub -cwd -V -b y -N \"S_${i}_${j}_${m}\" -l h_vmem=10G -l h_rt=24:00:00 -l h_stack=8M " \
#	"python2.7 mesa sedef/sedef search single fasta/hg19.fa chr$i chr$j $m"
```

