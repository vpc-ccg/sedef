# SEDEF: <u>Se</u>gmental <u>D</u>uplication <u>E</u>valuation <u>F</u>ramework

SEDEF is a quick tool to find all segmental duplications in the genome.

## Paper

SEDEF has been accepted at [ECCB 2018](http://eccb18.org). 
Preprint is [available here](https://arxiv.org/abs/1807.00205).

### Results

| 👨‍🎨 Human (hg19) | 🐭 Mouse (mm8) |
|-----|-----|
| [Final calls](http://cb.csail.mit.edu/cb/sedef/hg19.bed) | [Final calls](http://cb.csail.mit.edu/cb/sedef/mm8.bed) |

Paper experiments are outlined [in this Jupyter notebook](paper/experiments.ipynb).

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

> You need at least g++ 5.1.0 (C++14) to compile SEDEF. clang should work fine as well.

SEDEF requires Boost libraries in order to compile. In case you have non-standard Boost installation, you can still compile as follows:

```bash
CPATH={path_to_boost} make -j release
```

## How to run

The genome assembly must be soft-masked: all common and tandem repeats are converted to lower-case letters.
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
If you want to experiment with different parameters, run `sedef help` for parameter documentation.

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
		[ "$m" == "y" ] && rc="-r" || rc="";
		echo "sedef search $rc hg19.fa chr$i chr$j >out/${i}_${j}_${m}.bed 2>out/log/${i}_${j}_${m}.log"
		done; 
	fi
done
done | time parallel --will-cite -j 80 --eta

# Now make sure that all runs completed successfully 
grep Total out/log/*.log | wc -l
# You should see here 600 (or n(n+1) if you have n chromosomes in your file)

# Get the single-core running time
grep Wall out/log/*.log | tr -d '(' | awk '{s+=$4}END{print s}'

# Get the maximum meory usage as well
grep Memory out/log/*.log | awk '{if($3>m)m=$3}END{print m}'
```

Then use `sedef-align` to bucket the files for the optimal parallel alignment, and
afterwards run the whole alignment:
```bash
# First bucket the reads into 1000 bins
mkdir -p out/bins
mkdir -p out/log/bins
time sedef align bucket -n 1000 out out/bins

# Now run the alignment
for j in out/bins/bucket_???? ; do
	k=$(basename $j);
	echo "sedef align generate -k 11 hg19.fa $j >${j}.bed 2>out/log/bins/${k}.log"
done | time parallel --will-cite -j 80 --eta

# Make sure that all runs finished nicely
grep Finished out/log/bins/*.log | wc -l
# Should be number of bins (in our case, 1000)

# Get again the total running time
grep Wall out/log/bins/*.log | tr -d '(' | awk '{s+=$4}END{print s}'

# And the memory
grep Memory out/log/bins/*.log | awk '{if($3>m)m=$3}END{print m}'
```

Finally, run `sedef-stats` to produce the final output:
```bash
# Concatenate the files 
cat out/*.bed > out.bed # seed SDs
cat out/bins/bucket_???? > out.init.bed # potential SD regions
cat out/bins/*.bed | sort -k1,1V -k9,9r -k10,10r -k4,4V -k2,2n -k3,3n -k5,5n -k6,6n |\
	uniq > out.final.bed # final chains

# Now get the final calls
sedef stats generate hg19.fa out.final.bed |\
		sort -k1,1V -k9,9r -k10,10r -k4,4V -k2,2n -k3,3n -k5,5n -k6,6n |\
		uniq > out.hg19.bed
```

## Acknowledgements 

SEDEF uses [{fmt}](https://github.com/fmtlib/fmt), [argh](https://github.com/adishavit/argh) and the modified version of [Heng Li's ksw2](https://github.com/lh3/ksw2).

## Support

Questions, bugs? Open a GitHub issue or drop me an e-mail at `inumanag at mit dot edu`.

