#!/bin/bash

for fac_dens in 1 2 4 8 16 32 64 128 256 512; do
	echo "${fac_dens}"
	./bin/rosa_sd2_delta --input_file=web-4GB --add_factor_border --fac_dens=${fac_dens}
done
