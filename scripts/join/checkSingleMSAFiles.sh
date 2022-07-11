#!/bin/bash

a3mdir=$1
for f in $a3mdir/*.a3m
do
   fn=$(basename $f)
   name=${fn%.*}
   linecount=$(wc -l < $a3mdir/$name\.a3m)
   #lin=$(($linecount))
   #tt=$((lin))
   #$echo "$((tt))"
   #force even number of lines
   if [ $linecount = 2 ]; then
       echo "HERE!!! $linecount"
       echo "$name" >> one_msa.txt
   fi

done
