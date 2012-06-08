#!/bin/bash

bin_dir="../../../bin"
output_dir="../../../data/output/indexes"
tmp_dir="../../../data/tmp"
file="/devhome/sgog/rosa/paper/experiments/web-4GB"
file_basename=`basename ${file}`

threshold=4096

pattern_number=1000
pattern_len=20
result_file="query_performance_web-4GB.txt"
repetitions=5

rm -f ${result_file}

for i in `seq 1 ${repetitions}`; do
	for pattern_len in 4 10 20 40 100; do
		for occ in "1 2" "8 12" "75 125" "750 1250" "7500 12500"; do 
			min_occ=${occ% *}
			max_occ=${occ#* }
			echo "pattern_len=${pattern_len} min_occ=${min_occ} max_occ=${max_occ}"
			pattern_file=${output_dir}/${file_basename}.${pattern_len}.${pattern_number}.0.${min_occ}.${max_occ}.pattern
			echo "${pattern_file}"
			${bin_dir}/rosa_sd --input_file=${file} --threshold=${threshold} --output_dir=${output_dir} \
							   --benchmark --pattern_file=${pattern_file}  >> ${result_file}
		done
	done
done
