#!/bin/bash

root_dir=$(dirname $0)
download_dir=$root_dir/databases

if [ ! -d $download_dir ]; then
    mkdir -p $download_dir
fi
echo "Downloading BFD database"
bfd="$download_dir/bfd"
mkdir -p "$bfd"
bfd_filename="bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz"
wget -P "$bfd" "https://storage.googleapis.com/alphafold-databases/casp14_versions/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz"
if [ $?!=0 ]; then
    echo "Error!!! Unable to download BFD. Quitting."
    exit 11
fi
echo "BFD successfully downloaded. Now extracting the tarball ..."
tar --extract --verbose --file="$bfd/$bfd_filename" --directory="$bfd"
if [ $?!=0 ]; then
    echo "Error!!! Unable to decompress BFD. Quitting."
    exit 12
fi
echo "BFD extracted ... cleaning up."
rm "$bfd/$bfd_filename"

echo "Downloading Modeller database"
pir="$download_dir/pir"
mkdir -p "$pir"
pir_filename="pdball.pir.gz"
wget -P "$pir" "https://salilab.org/modeller/downloads/pdball.pir.gz"
if [ $?!=0 ]; then
    echo "Error!!! Unable to download Modeller PIR file. Quitting."
    exit 22
fi
echo "Modeller PIR successfully downloaded. Now extracting the tarball ..."
gunzip $pir/$pir_filename 
if [ $?!=0 ]; then
    echo "Error!!! Unable to decompress PIR file. Quitting."
    exit 23
fi

echo "Modeller PIR extracted ... creating hdf5 file."
conda activate dagmsa_env
if [ $?!=0 ]; then
    echo "Error!!! Unable to activate 'dagmsa_env'. Quitting."
    exit 24
fi
python $root_dir/scrtips/format_pir2hdf5.py $pir/pdball.pir $pir/
if [ $?!=0 ]; then
    echo "Error!!! Unable to create '$pir/pdball.hdf5' file. Some stuff may not run."
    #exit 24
fi
echo "Database downloads are complete!"