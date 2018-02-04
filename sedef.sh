#!/usr/bin/env bash
# 786

getopt --test > /dev/null
if [[ $? -ne 4 ]]; then
	echo "I’m sorry, `getopt --test` failed in this environment."
	exit 1
fi

PATH="${PATH}:"`pwd`

if ! command -v "samtools" >/dev/null 2>&1 ; then
	echo "Samtools not found in \$PATH"
	exit 1
fi

if ! command -v "parallel" >/dev/null 2>&1 ; then
	echo "GNU Parallel not found in \$PATH"
	exit 1
fi

if ! command -v "sedef" >/dev/null 2>&1 ; then
	echo "SEDEF not found in \$PATH (${PATH})"
	exit 1
fi

OPTIONS=hj:o:f
LONGOPTIONS=help,jobs,output,force
PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@")
if [[ $? -ne 0 ]]; then
	exit 2
fi
eval set -- "$PARSED"

output="sedef_out"
jobs=4
force="n"
while true; do
	case "$1" in
		-h|--help)
			show_help
			exit 0
			;;
		-f|--force)
			force="y"
			shift
			;;
		-o|--output)
			output="$2"
			shift 2
			;;
		 -j|--jobs)
			jobs="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			echo "Programming error"
			exit 3
			;;
	esac
done

if [[ $# -ne 1 ]]; then
	echo "$0: FASTA file is required."
	exit 1
fi
input="$1"

echo "SEDEF: FASTA=${input}; output=${output}; jobs=${jobs}; force=${force}"

if [ ! -f "${input}" ]; then
    echo "File ${input} not found!"
    exit 1
fi

if [ -e "${output}" ]; then
    echo -n "Output file name ${output} exists!"
    if [ "${force}" == "y" ] ; then
    	echo " Removing it."
    	rm -rf "${output}"
    else
    	echo " Please delete ${output} or run with -f/--force to proceed."
    	exit 1
    fi
fi
mkdir -p "${output}"


if [ ! -f "${input}.fai" ]; then
    echo "Indexing ${input}..."
    samtools faidx "${input}"
fi

mkdir -p "${output}/seeds"
mkdir -p "${output}/log/seeds"

echo "************************************************************************"
echo "Running SD seeding..."
for i in `cut -f1 "${input}.fai"`; do 
	for j in `cut -f1 "${input}.fai"`; do  
		SI=`awk '$1=="'$i'" {print $2}' "${input}.fai"`
		SJ=`awk '$1=="'$j'" {print $2}' "${input}.fai"` 
		if [ "$SI" -le "$SJ" ] ; then 
			for m in y n ; do
				echo "/usr/bin/time -f'TIMING: %e %M' sedef search single ${input} $i $j $m >${output}/seeds/${i}_${j}_${m}.bed 2>${output}/log/seeds/${i}_${j}_${m}.log"
			done
		fi
	done
done | tee "${output}/seeds.comm" | /usr/bin/time -f'Seeding time: %E' parallel --will-cite -j ${jobs} --eta

proc=`cat "${output}/seeds.comm" | wc -l`
echo "SD seeding done: done running ${proc} jobs!"

proc_ok=`grep Total ${output}/log/seeds/*.log | wc -l`
if [ "${proc}" != "${proc_ok}" ]; then
	echo "Error: launched ${proc} jobs but completed only ${proc_ok} jobs; exiting..."
	exit 2
fi

# Get the single-core running time
sc_time=`grep TIMING ${output}/log/seeds/*.log | awk '{s+=$2}END{print s}'`
sc_time_h=`echo "${sc_time} / 3600" | bc`
echo "Single-core running time: ${sc_time_h} hours (${sc_time} seconds)"

sc_mem=`grep TIMING ${output}/log/seeds/*.log | awk '{if($3>m)m=$3}END{print m}'`
sc_mem_k=`echo "${sc_mem} / 1024" | bc`
echo "Memory used: ${sc_mem_k} MB"


echo "************************************************************************"
echo "Running SD alignment..."

mkdir -p "${output}/align"
mkdir -p "${output}/log/align"
/usr/bin/time -f'Bucketing time: %E' sedef align bucket "${output}/seeds" "${output}/align" 1000 2>"${output}/log/bucket.log"

# Now run the alignment
for j in "${output}/align/bucket_"???? ; do
	k=$(basename $j);
	echo "/usr/bin/time -f'TIMING: %e %M' sedef align generate \"${input}\" $j 11 >${j}.aligned.bed 2>${output}/log/align/${k}.log"
done | tee "${output}/align.comm" | /usr/bin/time -f'Aligning time: %E' parallel --will-cite -j "${jobs}" --eta

proc=`cat "${output}/align.comm" | wc -l`
echo "SD alignment done: finished ${proc} jobs!"

proc_ok=`grep Finished ${output}/log/align/*.log | wc -l`
if [ "${proc}" != "${proc_ok}" ]; then
	echo "Error: launched ${proc} jobs but completed only ${proc_ok} jobs; exiting..."
	exit 2
fi

# Get the single-core running time
sc_time=`grep TIMING ${output}/log/align/*.log | awk '{s+=$2}END{print s}'`
sc_time_h=`echo "${sc_time} / 3600" | bc`
echo "Single-core running time: ${sc_time_h} hours (${sc_time} seconds)"

sc_mem=`grep TIMING ${output}/log/align/*.log | awk '{if($3>m)m=$3}END{print m}'`
sc_mem_k=`echo "${sc_mem} / 1024" | bc`
echo "Memory used: ${sc_mem_k} MB"


echo "************************************************************************"
echo "Running SD reporting..."

cat "${output}/seeds/"*.bed > "${output}/seeds.bed" # seed SDs
cat "${output}/align/bucket_"???? > "${output}/potentials.bed" # potential SD regions
cat "${output}/align/"*.aligned.bed > "${output}/aligned.bed"  # final chains

wc -l "${output}/"*.bed

# Now get the final calls
/usr/bin/time -f'Report time: %E' sedef stats generate "${input}" "${output}/aligned.bed" > "${output}/final.bed"

echo "************************************************************************"
echo "SEDEF done! Final SDs available in ${output}/final.bed"


