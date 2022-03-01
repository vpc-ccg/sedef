# 🔴 ⚠️ SEDEF has been deprecated. Please use [BISER](https://github.com/0xTCG/biser/) (SEDEF's successor) instead. ⚠️ 🔴


# SEDEF: <u>Se</u>gmental <u>D</u>uplication <u>E</u>valuation <u>F</u>ramework

SEDEF is a tool for quick detection of segmental duplications in a genome.

## Paper

SEDEF has been presented at [ECCB 2018](http://eccb18.org) (DOI [10.1093/bioinformatics/bty586](https://doi.org/10.1093/bioinformatics/bty586)). 
Preprint is [available here](https://arxiv.org/abs/1807.00205).
Get the final paper [here](https://academic.oup.com/bioinformatics/article/34/17/i706/5093240).

### Results

| 👨🎨 Human (hg38) | 👨‍🎨 Human (hg19) | 🐭 Mouse (mm8) |
|-----|-----|-----|
| [Final calls](http://alkanlab.org/share/sedef/hg38.bed) | [Final calls](http://cb.csail.mit.edu/cb/sedef/hg19.bed) | [Final calls](http://cb.csail.mit.edu/cb/sedef/mm8.bed) |

The experiment pipeline from the paper is described [in this Jupyter notebook](paper/experiments.ipynb).

## How to compile

Simple! Do this:
```bash
git clone https://github.com/vpc-ccg/sedef
cd sedef
make -j release
```

By default, SEDEF uses Intel C++ compiler. If you are using g++, build with:
```bash
make -j release CXX=g++
```

If you are using Clang on macOS, compile as 
```bash
brew install libomp
make -j release OPENMP="-Xpreprocessor -fopenmp" CXX=clang++
```

> You need at least g++ 5.1.0 (C++14) to compile SEDEF. Clang should work fine as well.

SEDEF requires Boost libraries in order to compile. In case you installed Boost in a non-standard directory, you can still compile as follows:
```bash
CPATH={path_to_boost} make -j release
```

## How to run

The genome assembly **must be soft-masked** (i.e. all common and tandem repeats should be converted to lower-case letters) and **indexed**.
Suppose that our genome is `hg19.fa` (we use UCSC hg19 genome with 24 standard chromosomes that does not contain patches (unGl) or random strains (chrXX_random)).

### Automatic transmission

Just go to `sedef` directory and run
```bash
./sedef.sh -o <output> -j <jobs> <genome> 
```

For example, to run hg19.fa on 80 cores type:
```bash
./sedef.sh -o sedef_hg19 -j 80 hg19.fa
```

You can add `-f` if `sedef_hg19` already exists (it will overwrite the existing content though). 
The final results will be located in `sedef_hg19/final.bed`.

Please note that `sedef.sh` depends on Samtools and GNU Parallel.
If you want to experiment with different parameters, run `sedef help` for parameter documentation.

### Output

Output will be located in `<out_dir>/final.bed`.

The fields of BEDPE file are as follows:

First 6 fields are standard BEDPE fields describing the coordinates of SD mates:
- `chr1`, `start1` and `end1`
- `chr2`, `start2` and `end2`

Other fields are (in the order of appearance):

| Field              | Description |
|--------------------|--------------------|
| `name`             | SD name |
| `score`            | Total alignment error  |
| `strand1`          | 1st SD mate strand |
| `strand2`          | 2nd SD mate strand |
| `max_len`          | Length of longer mate  |
| `aln_len`          | Alignment length (length with gaps) |
| `cigar`            | Empty string |
| `comment`          | Comment: currently shows mismatch base error (`m`) and gap base error (`g`) |
| `indel_a`          | Number of gap bases in the 1st mate |
| `indel_b`          | Number of gap bases in the 2nd mate     |
| `alnB`             | Aligned base count (matches and mismatches without gaps) |
| `matchB`           | Match base count  |
| `mismatchB`        | Mismatch base count |
| `transitionsB`     | Transition count (A <-> G and C <-> T) |
| `transversions`    | Transversion count (all mismatches that are not transitions)    |       
| `fracMatch`        | `matchB / alnB` |
| `fracMatchIndel`   | `matchB / aln_len`    |        
| `jck`              | Jaccard score: <img src="https://latex.codecogs.com/svg.latex?\frac{3}{4}\log\left(1-\frac{4}{3}w\right)" /> where `w = mismatchB / alnB` |
| `k2K`              | Kimura score: <img src="https://latex.codecogs.com/svg.latex?\frac{1}{2}\log\left(\frac{1}{1-2p-q}\right)+\frac{1}{4}\log\left(\frac{1}{1-2q}\right)" /> where `p = transitionsB / alnB` and `q = transitionsB / alnB` |
| `aln_gaps`         | Number of gaps in the alignment     |
| `uppercaseA`       | Number of non-masked (uppercase) bases in the 1st mate    |    
| `uppercaseB`       | Number of non-masked (uppercase) bases in the 2nd mate    |   
| `uppercaseMatches` | Non-masked match count |
| `aln_matches`      | Match base count |
| `aln_mismatches`   | Mismatch base count            |
| `aln_gaps`         | Number of gaps in the alignment      |
| `aln_gap_bases`    | Number of gap bases (`indel_a + indel_b`)    |     
| `cigar`            | CIGAR string of the SD mate alignment |
| `filter_score`     | `(aln_gaps + aln_mismatches) / aln_len` (should be ≥ 0.5) |

All errors are expressed as ratios (0.0--1.0) of the alignment length unless otherwise noted.

> *Warning*: as per WGAC, when calculating the similarity and error rates (fields `score`, `fracMatch`, 
`fracMatchIndel` and `filter_score`) SEDEF counts a gap as a *single* error 
(so the hypothetical alignment of `A-----GC` and `AT-----C` will have error 4 and NOT 8). This might lead
to SDs with rather large gap contents. For more filtering, consult `comment` field that provides the
percentage of match/mismatch and gap bases.

### Manual transmission

First make sure to index the genome:

```bash
samtools faidx hg19.fa
```

Then run the `sedef-search` in parallel (in this example, we will use GNU parallel) to get the initial seeds:
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

Then use `sedef-align` to bucket the files for the optimal parallel alignment. 
Afterwards, start the whole alignment:
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

## Acknowledgments 

SEDEF uses [{fmt}](https://github.com/fmtlib/fmt), [argh](https://github.com/adishavit/argh) and the modified version of [Heng Li's ksw2](https://github.com/lh3/ksw2).

## Support

Questions, bugs? Open a GitHub issue or drop me an e-mail at `inumanag at mit dot edu`.

