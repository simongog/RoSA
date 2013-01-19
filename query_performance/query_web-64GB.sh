#!/bin/bash

bin_dir="../bin"
output_dir=".."
pattern_dir="../pattern"
tmp_dir="../tmp"
file="../web-64GB"
file_basename=`basename ${file}`

threshold=4096

pattern_number=1000
pattern_len=20
result_file="query_performance_web-64GB.txt"
repetitions=5

if [[ $# -lt 1 ]] ; then
	echo "Please give the fac_dens as first command line argument "
	exit
fi 

fac_dens=${1}
os=`uname`

for i in `seq 1 ${repetitions}`; do
	for pattern_len in 40; do
		for occ in "8 12" "75 125"; do 
			if [[ "$os" == 'Darwin' ]]; then
				purge
			fi

			min_occ=${occ% *}
			max_occ=${occ#* }
			echo "pattern_len=${pattern_len} min_occ=${min_occ} max_occ=${max_occ}"
			pattern_file=${pattern_dir}/${file_basename}.${pattern_len}.${pattern_number}.0.${min_occ}.${max_occ}.pattern
			echo "${pattern_file}"
			${bin_dir}/rosa_sd2_delta --input_file=${file} --threshold=${threshold} --output_dir=${output_dir} \
							   --benchmark_ext --pattern_file=${pattern_file} --fac_dens=${fac_dens} >> ${result_file}
		done
	done
done
