#!/bin/bash
#this script will generate the hidden markov models of the combined MSA given no SS available.
#Usage:

#Copy the fastaA and fastaB to temp_dir as an a3m file
#use addss_v2020.pl to generate SS and add to the a3m file
#Read the SS from the ss.a3m file and write to top of fastaA_fastaB.a3m to make fastaA_fastaB.ss.a3m
fastafolder=/data/multicom4s_data/HETERO_STD/fastas
filename="5DQSA_5DQSD.a3m"
filename="${filename%.*}"
splits=$(echo $filename| tr "_" "\n")
outdir=/data/fqg7h/multicom4s_TMB/hmm_test
work_dir=$outdir/tmpdir$filename

#hhsuite3_dir=/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-3.0-beta.3-Linux #/home/multicom4s_tool/ZComplexMSA/bin/hhsuite-3.0-beta.3-Linux #/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-3.2.0-SSE2-Linux
#nr_db=/home/multicom4s_tool/alphafold_databases/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt  #### switch to BFD ###
#nr_iteration_num=4

#######
### SET #### Currently settling these ...
hhblits_dir=/home/fqg7h/multicom4s_TMB/tools_/hhblits/   #/home/fqg7h/multicom4s_TMB/tools/hhblits/   #### Change ####
multicom_tool_dir=/home/fqg7h/multicom4s_TMB/tools_/ ## Local ##### Change #####
multicom_database_dir=/home/fqg7h/multicom4s_TMB/databases ## Local ##### Change ####/data/farhan/multicom4s/tool

#### This is hhsuite 2.0 version 
hhsuite_dir=/data/farhan/multicom4s/tool/multicom4s_TMB/tools_/hhsuite-2.0.8-linux-x86_64/ #/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-2.0.16-linux-x86_64 #/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-3.0-beta.3-Linux  #/exports/store2/casp14/tools/hhsuite-2.0.8-linux-x86_64

#### Change the hhsuite DB ######
#### Currently settling these ...
hhsuitedb=/data/farhan/multicom4s/tool/multicom4s_TMB/tools_/update_db_hhsuite/hhsuitedb3 #/exports/store2/casp14/databases/update_db_hhsuite/hhsuitedb3 #####

##### Change this to one common directory
meta_dir=/data/farhan/multicom4s/tool/multicom4s_TMB/tools_/hhsuite
meta_common_dir=/exports/store2/casp14/tools/meta/


if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

if [ ! -d $work_dir ]; then
    mkdir -p $work_dir
fi


for names in $splits
do
    echo "$names"
    scp $fastafolder/$names\.fasta $work_dir/$names\.a3m
    echo "Adding Secondary structure information to a3m..."
    perl $hhblits_dir/scripts/addss_v2020.pl $multicom_tool_dir $multicom_database_dir $work_dir/$names.a3m $work_dir/$names.ss.a3m -a3m   
done