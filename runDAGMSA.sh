#!/bin/bash

#Usage: sh ./runDAGMSA.sh <fasta_list_file> <output_dir>

conda activate dagmsa_env
root_dir=$(readlink -f $0)
if [ ! -e $1 ]; then
	echo "Fasta list file not found!. Quitting"
	exit 1
fi

if [ ! -d $2 ]; then
	mkdir -p $2
fi

sh $root_dir/run_DAGMSA_pir.sh $1 $2
sh $root_dir/run_DAGMSA_hhb.sh $1 $2


