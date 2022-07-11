#!/bin/bash

# this script will take two fasta sequences of the dimers. Creates a concatenated MSA+SS and then generate ranked templates
source /home/multicom4s_tool/anaconda3/bin/activate
conda activate tbm_env

fasta_file_A=$1 #$(readlink - f $1)
fasta_file_B=$2 #$(readlink - f $1)
work_dir=$3 #$(readlink -f $2)

package_dir=$(dirname $(readlink -f $0))
echo "$package_dir"

#echo "$(readlink -f $fasta_file)"
fasta_file_A=$(readlink -f $fasta_file_A)
fasta_file_B=$(readlink -f $fasta_file_B)
work_dir=$(readlink -f $work_dir)
############## The following are the hard paths to the different tools. Change them accordingly. ####################

hhsuite3_dir=/home/multicom4s_tool/ZComplexMSA/bin/hhsuite-3.0-beta.3-Linux #/home/multicom4s_tool/ZComplexMSA/bin/hhsuite-3.0-beta.3-Linux #/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-3.2.0-SSE2-Linux
nr_db=/data/multicom4s_tool/alphafold_databases/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt  #### switch to BFD ###
nr_iteration_num=3
jack_dir=/home/multicom4s_tool/ZComplexMSA/bin/hmmer-3.1b2-linux-intel-x86_64

#######
### SET #### Currently settling these ...
hhblits_dir=$package_dir/tools_/hhblits/   #/home/fqg7h/multicom4s_TMB/tools/hhblits/   #### Change ####
multicom_tool_dir=$package_dir/tools_/ ## Local ##### Change #####
multicom_database_dir=$package_dir/databases ## Local ##### Change ####/data/farhan/multicom4s/tool
string_db=$package_dir/databases/string_db/strings_db_monomer_nr_1.00.txt
#### This is hhsuite 2.0 version 
hhsuite_dir=$package_dir/tools_/hhsuite-2.0.8-linux-x86_64/ #/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-2.0.16-linux-x86_64 #/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-3.0-beta.3-Linux  #/exports/store2/casp14/tools/hhsuite-2.0.8-linux-x86_64

#### Change the hhsuite DB ######
#### Currently settling these ...
hhsuitedb=$package_dir/tools_/update_db_hhsuite/hhsuitedb3 #/exports/store2/casp14/databases/update_db_hhsuite/hhsuitedb3 #####

##### Change this to one common directory
meta_dir=$package_dir/tools_/hhsuite
meta_common_dir=/exports/store2/casp14/tools/meta/

##### This is the .atom files of the homomers/heteromers 
#This needs to be fixed/set 
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
    echo "$fasta_file_A not found! Quitting!" #/data/multicom4s_tool/alphafold_databases/data/bfd
    exit 1
fi

if [ ! -e $fasta_file_B ]; then
    echo "$fasta_file_B not found! Quitting!"
    exit 1
fi


if [ ! -d $work_dir ]; then
    mkdir -p $work_dir
    chmod +777 $work_dir
fi

#work_dir=$(readlink -f $2)

#cp $fasta_file $work_dir/

cd $work_dir

### For MSA of Individual Fasta Files

scp $fasta_file_A $work_dir
scp $fasta_file_B $work_dir

scp $fasta_file_A $work_dir/$nameA.fasta.a3m
scp $fasta_file_B $work_dir/$nameB.fasta.a3m

############## Generate msa using hhblits ###########################
echo "Generate alignments from the bdf database using jackhmmer..."

echo "[Running] sh $package_dir/align_gen_hetero.sh $work_dir/$nameA.fasta $work_dir/$nameB.fasta $work_dir/ $string_db $jack_dir"

#####sh $package_dir/align_gen_hetero.sh $work_dir/$nameA.fasta $work_dir/$nameB.fasta $work_dir/ $string_db $jack_dir

#For First Fasta File
#$hhsuite3_dir/bin/hhblits -i $work_dir/$nameA.fasta -d $nr_db -oa3m $work_dir/$nameA.a3m -n $nr_iteration_num
if [ $? != 0 ]; then
    echo "Alignment generation using jackhmmer problem. Quitting !!!"
    exit 2
fi
#For Second Fasta File
#$hhsuite3_dir/bin/hhblits -i $work_dir/$nameB.fasta -d $nr_db -oa3m $work_dir/$nameB.a3m -n $nr_iteration_num
#if [ $? != 0 ]; then
#    echo "Alignment generation using hhblits problem for $fasta_file_B. Quitting !!!"
#    exit 2
#fi

echo "Alignments for $fasta_file_A and $fasta_file_B are successfully generated"

####### Concatenate the two MSAs ###################
echo "Concatenating the two alignment files ..."

#python $package_dir/scripts/joinAlignments4Homo_v2.py $work_dir/$nameA.a3m $work_dir/$nameB.a3m $work_dir

if [ $? != 0 ]; then
    echo "Unable to concatenate MSAs to generate paired MSA. Quitting !!!"
    exit 3
fi
echo "Paired MSA successfully added to $work_dir/$nameAB.a3m"


################ Add secondary structure into alignments #################
##### Modify the paths and add the tools in $hhblits_dir/scripts/addss_v2020.pl ###
echo "Adding Secondary structure information to a3m..."
#For First MSA
echo "[Running] perl $package_dir/perl_scripts/addss_v2020.pl $multicom_tool_dir $multicom_database_dir $work_dir/$nameA.fasta.a3m $work_dir/$nameA.fasta.ss.a3m -a3m"
#perl $package_dir/perl_scripts/addss_v2020.pl $multicom_tool_dir $multicom_database_dir $work_dir/$nameA.fasta.a3m $work_dir/$nameA.fasta.ss.a3m -a3m
perl $package_dir/tools_/hhblits/scripts/addss_v2020.pl $multicom_tool_dir $multicom_database_dir $work_dir/$nameA.fasta.a3m $work_dir/$nameA.fasta.ss.a3m -a3m

#exit
if [ $? != 0 ]; then
    echo "Unable to add Secondary structure information. Quitting !!!"
    exit 3
fi
echo "Secondary structure information successfully added to $work_dir/$nameA.fasta.ss.a3m"

#For second MSA
echo "perl $package_dir/tools_/hhblits/scripts/addss_v2020.pl $multicom_tool_dir $multicom_database_dir $work_dir/$nameB.fasta.a3m $work_dir/$nameB.fasta.ss.a3m -a3m"
perl $package_dir/tools_/hhblits/scripts/addss_v2020.pl $multicom_tool_dir $multicom_database_dir $work_dir/$nameB.fasta.a3m $work_dir/$nameB.fasta.ss.a3m -a3m

if [ $? != 0 ]; then
    echo "Unable to add Secondary structure information. Quitting !!!"
    exit 3
fi

echo "Secondary structure information successfully added to $work_dir/$nameB.fasta.ss.a3m"


exit

###### Concatenate the alignments into one file #########

#python $package_dir/scripts/combineSS_a3m_separatedimers.py $work_dir/$nameA.fasta.ss.a3m $work_dir/$nameB.fasta.ss.a3m $work_dir/

#python $package_dir/scripts/combineSS_a3m.py $work_dir/$nameA.fasta.ss.a3m $work_dir/$nameB.fasta.ss.a3m $work_dir/

#cat $work_dir/$nameAB.fasta.ss.a3m $work_dir/$nameAB.a3m > $work_dir/$nameAB.ss.a3m

#exit

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

#search hmm against the dimer database
echo "Searching $name.hmm against $hhsuitedb...\n";

#### Remove any previous search results ##########
if [ -f $work_dir/$name.hhr ]; then
    echo "Removing existing $work_dir/$name.hhr ..."
    #rm $work_dir/$name.hhr
fi
echo "[Running] $hhsuite3_dir/bin/hhsearch -i $work_dir/$nameAB.hmm -d $hhsuitedb -realign -mact 0"

$hhsuite3_dir/bin/hhsearch -i $work_dir/$nameAB.hmm -d $hhsuitedb -realign -mact 0

if [ ! -f $work_dir/$nameAB.hhr ]; then
    echo "$work_dir/$nameAB.hhr is not generated by hhsuite search. Try one more time using $nameAB.a3m ...\n"
    $hhsuite3_dir/bin/hhsearch -i $nameAB.ss.a3m -d $hhsuitedb -realign -mact 0
    if [ -f $work_dir/$nameAB.hhr ]; then
        mv $nameAB.hhr $nameAB.hhr
    fi
fi

if [ $? != 0 ]; then
    echo "Unable to perform hhsearch. Quitting !!!"
    exit 5
fi

echo "hhsearch successfully done"
##### End of Search #####

###### Now the rest #############

### 1. Rank the templates
echo "Generate ranking list..."
echo "[Running] $meta_dir/script/rank_templates.pl $nameAB.hhr $work_dir/$nameAB.rank"
$package_dir/perl_scripts/rank_templates.pl $nameAB.hhr $work_dir/$nameAB.rank

if [ $? != 0 ]; then
    echo "Unable to rank the templates. Quitting !!!"
    exit 5
fi

echo "Ranking successfully done"

### 3. Parse the blast output
echo "parse hhsearch output..."
$package_dir/perl_scripts/parse_hhsearch.pl $work_dir/$nameAB.hhr $work_dir/$nameAB.local
if [ $? != 0 ]; then
    echo "Unable to parse $nameAB.hhr file. Quitting !!!"
    exit 6
fi

######## Use $nameAB.local to generate top ten template .pir files
python $package_dir/scripts/convert2pir.py $work_dir/$nameAB.rank $work_dir/$nameAB.local $work_dir/$nameA.fasta $work_dir/$nameB.fasta $work_dir/


#TO DO ###### Facing issues with Modeller #######$work_dir/$work_dir/

### 4. Validating the local alignments
#$package_dir/perl_scripts/validate_local.pl $fasta_file.local $atom_dir $fasta_file.local

##### MULTICOM ######
##### Preprocess local alignments. Address domain issues.
#$meta_common_dir/script/local_global_align.pl $align_option $fasta_file.local $fasta_file $fasta_file
#if [ $? != 0 ]; then
#    echo "Unable to process local alignments. Quitting !!!"
#    exit 7
#fi

#### Generate Model #####
#$meta_common_dir/script/local2model.pl $option_file $fasta_file $fasta_file $work_dir




