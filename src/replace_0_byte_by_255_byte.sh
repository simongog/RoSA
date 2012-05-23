#!/bin/bash
# Author: Simon Gog (simon.gog@unimelb.edu.au)

if [ $# -le 1 ];
then
	echo "usage: $0 original_file destination_file"
	echo " Replaces all 0-bytes in original_file by 255-bytes "
	echo " and writes the result into destination_file"
else
	tr '\000' '\377' < $1 > $2
fi

