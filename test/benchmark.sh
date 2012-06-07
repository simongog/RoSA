#!/bin/bash

for len in 5 10 20; do 
	PATTERN_NAME=english50MB.pattern
	../bin/rosa_sd --input_file=$1 --threshold=40000 --generate_patterns --pattern_len=${len} --pattern_number=10000 --output_dir=../data/output/indexes --pattern_file=${PATTERN_NAME}
	../bin/rosa_sd --input_file=$1 --threshold=40000 --benchmark_ext --pattern_file=${PATTERN_NAME} --output_dir=../data/output/indexes/ 
done

