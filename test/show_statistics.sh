#!/bin/bash

../bin/rosa_sd --input_file=$1 --threshold=4000 --output_statistics  --output_dir=../data/output/indexes
../bin/rosa_sd --input_file=$1 --threshold=40000 --output_statistics  --output_dir=../data/output/indexes
