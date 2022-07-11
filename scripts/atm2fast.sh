#!/bin/bash

infolder=$(readlink -f $1)
outfolder=$(readlink -f $2)

if [ ! -d $outfolder ]; then
   mkdir -p $outfolder
fi

for f in $infolder/*.atm; do
   fn=$(basename $f)
   name=${fn%.*}
   echo ">$name" > $outfolder/$name\.fasta
   /data/farhan/multicom4s/tool/multicom4s_TMB/scripts/pdb2seq.pl $f >> $outfolder/$name\.fasta
done
