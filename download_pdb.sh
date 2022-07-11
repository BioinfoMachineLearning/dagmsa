#!/bin/bash

list_file=$(readlink -f $1)
download_dir=$2
download_dir=$(readlink -f $download_dir)

if [ ! -d $download_dir ]; then
	mkdir -p $download_dir
fi

cd $download_dir

extension=$3

while IFS= read -r line
do
	echo "$line"
	echo "wget https://files.rcsb.org/download/$line.$extension"
	wget https://files.rcsb.org/download/$line.$extension
	mv $download_dir/$line.$extension $download_dir/$line\.pdb
	
	#break
done < "$list_file"
