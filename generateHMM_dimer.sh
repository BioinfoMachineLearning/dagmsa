#!/bin/bash
# this script will run the makeHHB.sh on the two fasta sequences of the dimers. Creates a concatenated MSA+SS
fasta_file_A=$1 #$(readlink - f $1)
fasta_file_B=$2 #$(readlink - f $1)
work_dir=$3 #$(readlink -f $2)

#echo "$(readlink -f $fasta_file)"
fasta_file_A=$(readlink -f $fasta_file_A)
fasta_file_B=$(readlink -f $fasta_file_B)
work_dir=$(readlink -f $work_dir)
############## The following are the hard paths to the different tools. Change them accordingly. ####################

hhsuite3_dir=/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-3.0-beta.3-Linux #/home/multicom4s_tool/ZComplexMSA/bin/hhsuite-3.0-beta.3-Linux #/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-3.2.0-SSE2-Linux
nr_db=/home/multicom4s_tool/alphafold_databases/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt  #### switch to BFD ###
nr_iteration_num=4

#######
### SET #### Currently settling these ...
hhblits_dir=/data/farhan/multicom4s/tool/multicom4s_TMB/tools_/hhblits/   #/home/fqg7h/multicom4s_TMB/tools/hhblits/   #### Change ####
multicom_tool_dir=/data/farhan/multicom4s/tool/multicom4s_TMB/tools_/ ## Local ##### Change #####
multicom_database_dir=/data/farhan/multicom4s/tool/multicom4s_TMB/databases ## Local ##### Change ####/data/farhan/multicom4s/tool

#### This is hhsuite 2.0 version 
hhsuite_dir=/data/farhan/multicom4s/tool/multicom4s_TMB/tools_/hhsuite-2.0.8-linux-x86_64/ #/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-2.0.16-linux-x86_64 #/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-3.0-beta.3-Linux  #/exports/store2/casp14/tools/hhsuite-2.0.8-linux-x86_64

#### Change the hhsuite DB ######
#### Currently settling these ...
hhsuitedb=/data/farhan/multicom4s/tool/multicom4s_TMB/tools_/update_db_hhsuite/hhsuitedb3 #/exports/store2/casp14/databases/update_db_hhsuite/hhsuitedb3 #####

##### Change this to one common directory
meta_dir=/data/farhan/multicom4s/tool/multicom4s_TMB/tools_/hhsuite
meta_common_dir=/exports/store2/casp14/tools/meta/

##### This is the .atom files of the heteromers
atom_dir=/exports/store2/casp14/databases/RCSB_PDB/atom/
align_option=$meta_dir/align_option  # Change the info at the align option file 
option_file=$same #Same stuff as here
prosys_dir=/exports/store2/casp14/tools/prosys/
script_dir=$prosys/scriptsacct

####################### End of Hard paths ##################################

fnA=$(basename $fasta_file_A)
nameA=${fnA%.*}
fnB=$(basename $fasta_file_B)
nameB=${fnB%.*}

nameAB=$nameA\_$nameB

echo "FASTA_A $fasta_file_A"
echo "NAME_A $nameA"
echo "FASTA_B $fasta_file_B"
echo "NAME_A $nameB"
echo "WORK_DIR $work_dir"
echo "NAME_AB $nameAB"
#exit 10

if [ ! -e $fasta_file_A ]; then
    echo "$fasta_file_A not found! Quitting!"
    exit 1
fi

if [ ! -e $fasta_file_B ]; then
    echo "$fasta_file_B not found! Quitting!"
    exit 1
fi


if [ ! -d $work_dir ]; then
    mkdir -p $work_dir
    chmod +775 $work_dir
fi

#work_dir=$(readlink -f $2)

#cp $fasta_file $work_dir/

cd $work_dir

### For MSA of Individual Fasta Files

scp $fasta_file_A $work_dir
scp $fasta_file_B $work_dir

############## Generate msa using hhblits ###########################
#echo "Generate alignments from the bdf database using hhblits..."
#For First Fasta File
#$hhsuite3_dir/bin/hhblits -i $work_dir/$nameA.fasta -d $nr_db -oa3m $work_dir/$nameA.a3m -n $nr_iteration_num
#if [ $? != 0 ]; then
#    echo "Alignment generation using hhblits problem for $fasta_file_A. Quitting !!!"
#    exit 2
#fi
#For Second Fasta File
#$hhsuite3_dir/bin/hhblits -i $work_dir/$nameB.fasta -d $nr_db -oa3m $work_dir/$nameB.a3m -n $nr_iteration_num
#if [ $? != 0 ]; then
#    echo "Alignment generation using hhblits problem for $fasta_file_B. Quitting !!!"
#    exit 2
#fi

#echo "Alignments for $fasta_file_A and $fasta_file_B are successfully generated"

################ Add secondary structure into alignments #################
##### Modify the paths and add the tools in $hhblits_dir/scripts/addss_v2020.pl ###
echo "Adding Secondary structure information to a3m..."
#For First MSA

##perl $hhblits_dir/scripts/addss_v2020.pl $multicom_tool_dir $multicom_database_dir $work_dir/$nameA.a3m $work_dir/$nameA.ss.a3m -a3m

if [ $? != 0 ]; then/data/farhan/SoftwareTools/DeepComplex/scripts/database_generating_scripts/hetero_std/pdb
if [ $? != 0 ]; then
    echo "Unable to add Secondary structure information. Quitting !!!"
    exit 3
fi

echo "Secondary structure information successfully added to $work_dir/$nameB.ss.a3m"

###### Concatenate the alignments into one file #########

##python /data/farhan/multicom4s/tool/multicom4s_TMB/scripts/joinAlignments_v1.py $work_dir/$nameA.a3m $work_dir/$nameB.a3m $work_dir
##python /data/farhan/multicom4s/tool/multicom4s_TMB/scripts/joinAlignments_ss_v1.py $work_dir/$nameA.ss.a3m $work_dir/$nameB.ss.a3m $work_dir



#### This is hhsuite 2.0 version 
#export HHLIB="$hhsuite_dir/lib/hh"


export HHLIB="$hhsuite3_dir"
##### Create the hidden markov model #############
echo "Creating Hidden Markov Models..." 
echo "[Running] $hhsuite3_dir/bin/hhmake -i $nameAB.ss.a3m -o $nameAB.hmm"
$hhsuite3_dir/bin/hhmake -i $nameAB.ss.a3m -o $nameAB.hmm
if [ $? != 0 ]; then
    echo "Unable to create HMM. Quitting !!!"
    #exit 4
fi

echo "HMM successfully created $workd_dir/$nameAB.hmm"

#search hhm against the database
echo "Searching $nameAB.hmm against $hhsuitedb...\n";

