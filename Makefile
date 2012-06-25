INCLUDE_PATH=${HOME}/include
LIB_PATH=${HOME}/lib
SRC_DIR=src
BIN_DIR=bin
CFLAGS=-O3 -funroll-loops -msse4.2 -DMEM_INFO -DOUTPUT_STATS -DMEM_INFO -DWRITE_R_OUTPUT -DNDEBUG # -DSDSL_DEBUG_ALGORITHMS_FOR_COMPRESSED_SUFFIX_ARRAYS

all: ${BIN_DIR}/rosa_helping_structures.o \
	 ${BIN_DIR}/rosa_helping_functions.o \
	 ${BIN_DIR}/pattern_file.o \
	 ${BIN_DIR}/rosa_sd \
	 ${BIN_DIR}/rosa_sd_load_only\
	 ${BIN_DIR}/rosa_sd_create_only\
	 ${BIN_DIR}/rosa_sd_search_block_only
#	 ${BIN_DIR}/rosa_bv \
#	 ${BIN_DIR}/rosa_rrr \
#	 ${BIN_DIR}/rosa_rrr63 \
#	 ${BIN_DIR}/rosa_gap \
#	 ${BIN_DIR}/fm_huff \
#	 ${BIN_DIR}/fm_huff_rrr \
#	 ${BIN_DIR}/fm_huff_rrr127 \
#	 ${BIN_DIR}/fm_rlmn

${BIN_DIR}/rosa_helping_structures.o:
	g++ ${CFLAGS} -I${INCLUDE_PATH} -c ${SRC_DIR}/rosa_helping_structures.cpp -o ${BIN_DIR}/rosa_helping_structures.o

${BIN_DIR}/rosa_helping_functions.o:
	g++ ${CFLAGS} -I${INCLUDE_PATH} -c ${SRC_DIR}/rosa_helping_functions.cpp -o ${BIN_DIR}/rosa_helping_functions.o

${BIN_DIR}/pattern_file.o:
	g++ ${CFLAGS} -I${INCLUDE_PATH} -c ${SRC_DIR}/pattern_file.cpp -o ${BIN_DIR}/pattern_file.o

${BIN_DIR}/rosa_sd: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} \
	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/rosa_sd \
		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
		-lsdsl -ldivsufsort -ldivsufsort64

${BIN_DIR}/rosa_sd_load_only: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DBENCHMARK_LOAD_ONLY \
	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/rosa_sd_load_only \
		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
		-lsdsl -ldivsufsort -ldivsufsort64

${BIN_DIR}/rosa_sd_create_only: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DBENCHMARK_CREATE_ONLY \
	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/rosa_sd_create_only \
		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
		-lsdsl -ldivsufsort -ldivsufsort64

${BIN_DIR}/rosa_sd_search_block_only: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DBENCHMARK_SEARCH_BLOCK_ONLY \
	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/rosa_sd_search_block_only \
		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
		-lsdsl -ldivsufsort -ldivsufsort64

#${BIN_DIR}/rosa_bv: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
#	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DROSA_BV \
#	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/rosa_bv \
#		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
#		-lsdsl -ldivsufsort -ldivsufsort64
#
#${BIN_DIR}/rosa_rrr: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
#	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DROSA_RRR \
#	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/rosa_rrr \
#		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
#		-lsdsl -ldivsufsort -ldivsufsort64
#
#${BIN_DIR}/rosa_gap: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
#	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DROSA_GAP \
#	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/rosa_gap \
#		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
#		-lsdsl -ldivsufsort -ldivsufsort64
#
#${BIN_DIR}/rosa_rrr63: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
#	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DROSA_RRR_VAR \
#	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/rosa_rrr63 \
#		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
#		-lsdsl -ldivsufsort -ldivsufsort64
#
#${BIN_DIR}/fm_huff: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
#	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DFM_HUFF \
#	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/fm_huff \
#		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
#		-lsdsl -ldivsufsort -ldivsufsort64
#
#${BIN_DIR}/fm_huff_rrr: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
#	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DFM_HUFF_RRR \
#	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/fm_huff_rrr \
#		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
#		-lsdsl -ldivsufsort -ldivsufsort64
#
#${BIN_DIR}/fm_huff_rrr127: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
#	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DFM_HUFF_RRR127 \
#	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/fm_huff_rrr127 \
#		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
#		-lsdsl -ldivsufsort -ldivsufsort64
#
#${BIN_DIR}/fm_rlmn: ${BIN_DIR}/pattern_file.o ${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o
#	g++ ${CFLAGS} -I${INCLUDE_PATH} -L${LIB_PATH} -DFM_RL \
#	    ${SRC_DIR}/rosa_main.cpp -o ${BIN_DIR}/fm_rlmn \
#		${BIN_DIR}/rosa_helping_functions.o ${BIN_DIR}/rosa_helping_structures.o ${BIN_DIR}/pattern_file.o \
#		-lsdsl -ldivsufsort -ldivsufsort64
#
clean:
	rm -f ${BIN_DIR}/*
