#!/bin/bash
#Usage: sh install_dagmsa.sh --db 
#This will install dagmsa along with databases
root_dir=$(readlink -f $0)
echo "Attempting to create conda environment 'dagmsa_env' ..."
conda create --name dagmsa_env python==3.8

if [ $? != 0 ]; then
   echo "Enable to create conda environment. Quitting!!!"
   exit 1
fi
echo "conda environment 'dagmsa_env' successfully created."


conda activate dagmsa_env

if [ $? != 0 ]; then
   echo "Enable to activate conda environment 'dagmsa_env'. Quitting!!!"
   exit 2
fi

echo "Installing Biopython ..."

conda install -c conda-forge biopython

if [ $? != 0 ]; then
   echo "Enable to install 'Biopython'. Some stuff may not run!!!"
   #exit 2
fi

echo "Biopython successfully installed ..."

echo "Installing HH-suite ..."

conda install -c bioconda hhsuite==3.3.0
if [ $? != 0 ]; then
   echo "Enable to install 'HH-suite'. Some stuff may not run!!!"
   #exit 2
fi

echo "HH-suite successfully installed"

echo "Installing Modeller ..."
#export KEY_MODELLER="MODELIRANJE"
conda install -c salilab modeller

if [ $? != 0 ]; then
   echo "Enable to install 'Modeller'. Some stuff may not run!!!"
   #exit 2
fi

if [ $1=="--db" ]; then
    echo "Now installing databases..."
    sh $0/download_db.sh 
fi

readlink -f $root_dir/example/1A1M/* > $root_dir/example/1A1M_list.txt

echo "Installation complete!"
