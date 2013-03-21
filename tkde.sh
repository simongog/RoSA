#!/bin/bash

for i in `seq 1 6`; do
for patsuf in "4.1000.0.7500.12500.pattern" "10.1000.0.750.1250.pattern" "20.1000.0.75.125.pattern" "40.1000.0.8.12.pattern" "100.1000.0.1.1.pattern"; do
	./bin/fm64_huff_rrr --fac_dens=16 --input_file=web-64GB --benchmark_ext --pattern_file=pattern/web-64GB.${patsuf}
done
done
