#!/bin/bash 
for fac_dens in 1 2 4 8 16 32 64 128 256 512 1024 0; do
	echo "fac_dens=${fac_dens}"
	./bin/rosa_sd2_delta --input_file=web-4GB --verbose --tmp_file_dir=tmp/ --fac_dens=${fac_dens} --generate_index
done
