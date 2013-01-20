#!/bin/bash

bin_dir="../bin"
output_dir=".."
pattern_dir="../pattern"
tmp_dir="../tmp"
file="../web-64GB"
file_basename=`basename ${file}`

threshold=4096

pattern_number=1000
result_file="query_performance_loc_web-64GB.txt"
repetitions=5

# rm -f ${result_file}

os=`uname`
#if [[ "$os" == 'Darwin' ]]; then
#	echo "Setup Mac OS for experiments"
#	./mac_setup.sh
#fi

#for fac_dens in 0 1 2 4 8 16 32 64 128 256 512 1024; do
for fac_dens in 0; do
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
					echo "${pattern_file} "
					echo "# phase = ${suffix}" >> ${result_file}
					${bin_dir}/rosa_sd2_delta --input_file=${file} --threshold=${threshold} --output_dir="${output_dir}/" \
								  --benchmark_loc --benchmark_ext --pattern_file=${pattern_file} --fac_dens=${fac_dens} >> ${result_file}
				done
		done
	done
done

#if [[ "$os" == 'Darwin' ]]; then
#	echo "Activate services on Mac OS again"
#	./mac_teardown.sh
#fi
