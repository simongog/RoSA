#!/bin/bash

../bin/rosa_sd --input_file=$1 --threshold=4096 --output_statistics  --output_dir=../data/output/indexes
../bin/rosa_sd --input_file=$1 --threshold=16384 --output_statistics  --output_dir=../data/output/indexes
