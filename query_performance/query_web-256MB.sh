#!/bin/bash

bin_dir="../bin"
output_dir=".."
tmp_dir="../tmp"
file="../web-256MB"
file_basename=`basename ${file}`

threshold=4096

pattern_number=1000
pattern_len=20
result_file="query_performance_web-256MB.txt"
repetitions=5
fac_dens=1

rm -f ${result_file}

for i in `seq 1 ${repetitions}`; do
	for pattern_len in 4 10 20 40 100; do
		for occ in "1 1" "8 12" "75 125" "750 1250" "7500 12500"; do 
			min_occ=${occ% *}
			max_occ=${occ#* }
			echo "pattern_len=${pattern_len} min_occ=${min_occ} max_occ=${max_occ}"
			pattern_file=${output_dir}/${file_basename}.${pattern_len}.${pattern_number}.0.${min_occ}.${max_occ}.pattern
			echo "${pattern_file}"
			${bin_dir}/rosa_sd2_delta --input_file=${file} --threshold=${threshold} --output_dir=${output_dir} \
							   --benchmark --pattern_file=${pattern_file} --fac_dens=${fac_dens} >> ${result_file}
		done
	done
done
