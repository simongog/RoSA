#!/bin/bash

threshold=3
fac_dens=1
rosa_exec=../bin/rosa_sd2_delta
tmp_dir=tmp

if [ $# -lt 1 ]; then
	echo "Usage: ${0} file [threshold] [fac_dens]" 
	echo "  file     : File containing the example string"
	echo "  threshold: Block threshold; default=${threshold}"
	echo "  fac_dens : Sampling parameter for factorization pointer; default=${fac_dens}"
	exit 1
fi

basename=`basename ${1}`

${rosa_exec} --input_file=${1} --threshold=${threshold} --fac_dens=${fac_dens} --generate_index --tmp_file_dir=${tmp_dir}/ --output_dir=${tmp_dir}/
${rosa_exec} --input_file=${1} --threshold=${threshold} --fac_dens=${fac_dens} --output_tikz --tmp_file_dir=${tmp_dir}/ --output_dir=${tmp_dir}/

tmp_tex_file=${basename}.${threshold}.tex
tmp_pdf_file=${basename}.${threshold}.pdf

sed "s|XXXXXXXXXX|${basename}.${threshold}|g" template/x.tex > ${tmp_dir}/${tmp_tex_file}

cd ${tmp_dir}
pdflatex ${tmp_tex_file}
pdflatex ${tmp_tex_file}
mv ${tmp_pdf_file} ..
cd ..

 rm ${tmp_dir}/*
