#!/bin/bash

bin_dir="../bin"
output_dir=".."
usb_dir="/Volumes/USB"
pattern_dir="../pattern"
tmp_dir="../tmp"
file="${usb_dir}/web-4GB"
file_basename=`basename ${file}`

threshold=4096

pattern_number=1000
pattern_len=20
result_file="query_performance_usb_web-4GB.txt"
repetitions=5

rm -f ${result_file}

os=`uname`
if [[ "$os" == 'Darwin' ]]; then
	echo "Setup Mac OS for experiments"
	./mac_setup.sh
fi

#for fac_dens in 0 1 2 4 8 16 32 64 128 256 512 1024; do
#for fac_dens in 0 1 4 16 64 256 1024; do
for fac_dens in 0; do
#	cp ${output_dir}/web-4GB.${threshold}.${fac_dens}.2.1.* ${usb_dir}
	for i in `seq 1 ${repetitions}`; do
		for pattern_len in 4 10 20 40 100; do
			for occ in "1 1" "8 12" "75 125" "750 1250" "7500 12500"; do 
			if [[ "$os" == 'Darwin' ]]; then
				purge
			fi
				min_occ=${occ% *}
				max_occ=${occ#* }
				echo "pattern_len=${pattern_len} min_occ=${min_occ} max_occ=${max_occ}"
				pattern_file=${pattern_dir}/${file_basename}.${pattern_len}.${pattern_number}.0.${min_occ}.${max_occ}.pattern
				echo "${pattern_file}"
				${bin_dir}/rosa_sd2_delta --input_file=${file} --threshold=${threshold} --output_dir="${usb_dir}" \
								   --benchmark --pattern_file=${pattern_file} --fac_dens=${fac_dens} >> ${result_file}
			done
		done
	done
#	rm ${usb_dir}/web-4GB.${threshold}.${fac_dens}.2.1.*
done

if [[ "$os" == 'Darwin' ]]; then
	echo "Activate services on Mac OS again"
	./mac_teardown.sh
fi
