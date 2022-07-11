#!/bin/bash

#Installation bash script for the external tools. This script already assumes packages are installed in the tools folder. It sets the required paths.
#NOTE: The paths for the database still needs to be set. It can be done through here or by separately updating the paths.txt file.
#NOTE: Only update the paths.txt file after running this code.
#1. gcc compiler. Assumes it is already installed. If not please install gcc version 4 or above separately before running this program. For clusters
#	type the following commands:
#	$ module avail gcc
#	chose the version you want to install. For example if version is gcc/gcc-5.4.0,
#	$ module load gcc/gcc-5.4.0
#2. cmake
#3. MPI 
#4. hh-suite
#5. MMSeqs
#6. JackHMMER
# Help or Usage message:

help_msg="""
usage: $0 
"""
packagedir=$(pwd)
mypath=$PATH

checkStatus()
{
	#msg = $1
	if [ ! $?==0 ]; then
		echo "Error code: $? $1"
		#rollBack $2
		#exit 1
	fi
}

echo "Loading module GCC-5.4.0..."
module load gcc/gcc-5.4.0

echo "Installing updated paths for the tools in the folder $packagedir/tools"

cd $packagedir/tools

cd hhsuite*
checkStatus "the directory for hhsuite was not found. Please check if tool is installed in $packagedir/tools folder"
HHLIB=$(pwd)

pathHhblits="$HHLIB/bin"
pathHhbdata="$HHLIB/data/"
pathHhbscripts="$HHLIB/scripts"

cd $packagedir/tools/hmmer*
checkStatus "the directory for hmmer was not found. Please check if tool is installed in $packagedir/tools folder"
pathJack=$(pwd)
pathJackbin=$(pwd)/src
patheasel=$(pwd)/easel/miniapps

#database paths please change accordingly
uniclustpath=/storage/htc/bdm/tools/databases/uniclust30_2017_10/
uniclustdb=uniclust30_2017_10
unirefpath=/storage/htc/bdm/tools/databases/uniref100_04_2018/
unirefdb=uniref100.fasta
ebipath=/storage/htc/bdm/tools/databases/unirefEBI/
ebi=ebimgnify.fasta
ebi_uniref_db=ebi_uniref100

cd $packagedir/tools/mmseqs*
checkStatus "the directory for mmseqs was not found. Please check if tool is installed in $packagedir/tools folder"
pathmmseqs_export=$(pwd)/bin
pathmmseqs=$(pwd)/bin/mmseqs

cd $packagedir/tools
pathmafft=$(pwd)/MAFFT/mafft-linux64/
mafft=mafft.bat

#hhblits coverage and threads. Change accordingly
cov=50
threads=8

cd $packagedir/tools/coneva
checkStatus "the directory for coneva was not found. Please check if tool is installed in $packagedir/tools folder"
pathconeva=$(pwd)
if [ ! -f $pathconeva/coneva.pl ]; then
	echo "$pathconeva/coneva.pl not found. Evaluation scripts will show errors"
fi

cd $packagedir/tools/CCMpred
checkStatus "the directory for CCMpred was not found. Please check if tool is installed in $packagedir/tools folder"
pathccmpred=$(pwd)/bin
if [ ! -f $pathccmpred/ccmpred ]; then
	echo "$pathccmpred/ccmpred not found. Evaluation scripts will show errors"
fi


cd $packagedir/tools/cmake*
checkStatus "the directory for cmake was not found. Please check if tool is installed in $packagedir/tools folder. Will not cause any execution problems. New installation might be a problem"
pathcmake=$(pwd)/bin
if [ ! -f $pathcmake/cmake ]; then
	echo "$pathcmake/cmake not found. Execution will not cause errors. New installation might be a problem"
fi



cd $packagedir/tools/mpi
checkStatus "the directory for openmpi was not found. Please check if tool is installed in $packagedir/tools folder. Will cause problems when generating alignments using jackhmmer."
pathmpi=$(pwd)/bin
pathmpi_lib=$(pwd)/lib

cd $packagedir
echo "#Change the paths accordingly" > temp.txt
echo "pathHhblits=$pathHhblits" >> temp.txt
echo "pathHhbdata=$pathHhbdata" >> temp.txt
echo "pathJack=$pathJack" >> temp.txt
echo "uniclustpath=$uniclustpath" >> temp.txt
echo "uniclustdb=$uniclustdb" >> temp.txt
echo "unirefpath=$unirefpath" >> temp.txt
echo "unirefdb=$unirefdb" >> temp.txt
echo "ebipath=$ebipath" >> temp.txt
echo "ebi=$ebi" >> temp.txt
echo "ebi_uniref_db=$ebi_uniref_db" >> temp.txt
echo "pathmmseqs=$pathmmseqs" >> temp.txt
echo "pathmafft=$pathmafft" >> temp.txt
echo "mafft=$mafft" >> temp.txt
echo "cov=$cov" >> temp.txt
echo "threads=$threads" >> temp.txt
echo "exports=export HHLIB=$HHLIB && export PATH=\$PATH:$pathHhblits:$pathHhbscripts:$pathJackbin:$patheasel:$pathconeva:$pathcmake:$pathmpi:$pathccmpred:$pathmmseqs_export && export LD_PATH_LIBRARY=\$LD_PATH_LIBRARY:$pathmpi_lib && " >> temp.txt

rm -f $packagedir/paths.txt

mv $packagedir/temp.txt $packagedir/paths.txt


export PATH=$PATH:$pathcmake:$pathmpi
export LD_PATH_LIBRARY=$LD_PATH_LIBRARY:$pathmpi_lib
cd $packagedir/tools/hhsuite*
INSTALL_BASE_DIR=$(pwd)
rm -r build
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX="$INSTALL_BASE_DIR" ..
make install
echo "All tools paths have been updated."
echo "Creating new python virtual environment..."
echo "Checking python version"
python $packagedir/scripts/checkPythonVersion.py
#if [ ! $?==0 ]; then
#	echo "Attempting to switch to python version 3.6"
#        sh $pckagedir/change_to_python3.sh
#	if [ ! $?==0 ]; then
#		echo "Python version change to 3.6 failed! Please install python version 3.6 or above to run this software. Quitting!"
#		exit 1
#	fi
#fi

if [ ! -d $packagedir/bin ]; then
	mkdir $packagedir/bin
fi
#cd $packagedir/bin

if [ -d $packagedir/bin/penv3.6 ]; then
	echo "Removing existing python virtual environment..."
	rm -rf $packagedir/bin/penv3.6
fi

module load python/python-3.6.6-tk
python3 -m venv $packagedir/bin/penv3.6
checkStatus "Error when creating python virtual environment. Please create it manually."

echo "Python virtual environment 'penv3.6' successfully created in the $packagedir/bin folder"

source $packagedir/bin/penv3.6/bin/activate
pip install numpy
checkStatus "Could not install numpy. Please install in manually."
deactivate

echo $packagedir > $packagedir/packagedir.txt
echo "Installation of Sensalign_v1 complete!"
#echo "Now run 'sh generate_alignment_4_sequence.sh' or 'sh generate_alignment_4_dir.sh'"
echo "Now run 'sh generate_alignment.sh <fastafile or fastadirectory > <outputfolder>'"

#mpirun
#mmseqs
