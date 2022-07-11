#!/bin/bash
# this script will generate hhblits database from list of hmm files
hhblits_dir=/data/farhan/SoftwareTools/DeepComplexToolBox/tools/hhsuite-3.0-beta.3-Source/bin #$1
hhm_dir=$1 #or msa_dir
outdir=$2

# concat all the a3m files into one file

if [ ! -d $outdir ]; then
    mkdir $outdir
fi

$hhblits_dir/ffindex_build -s $outdir/hhsuite3_db.ff{data,index}

$hhblits_dir/ffindex_apply $outdir/hhsuite3_db.ff{data,index} \ -i $outdir/hhsuite3_db_a3m.ffindex -d $outdir/hhsuite3_db_a3m.ffdata \ -- $hhblits_dir/hhconsensus -M 50 -maxres 65535 n-i stdin -oa3m -v 0
rm $outdir/hhsuite3_db.ff{data,index}
#$hhblits_dir/hhconcensus 

$hhblits_dir/ffindex_apply $outdir/hhsuite3_db_a3m.ff{ffdata,index} \\ -i $outdir/hhsuite3_db_hmm.ffindex -d hhsuite3_db_hmm.ffdata -- hhmake -i stdin -o stdout -v 0

$hhblits_dir/cstranslate -f -x 0.3 -c 4 I a3m -i $outdir/hhsuite3_db_a3m -o hhsuite3_db_cs219

sort -k3 -n -r $outdir/hhsuite3_db_cs219.ffindex | cut -f1 > sorting.dat

$hhblits_dir/ffindex_order sorting.dat $outdir/hhsuite3_db_hhm.ff{data,index} $outdir/hhsuite3_db_hhm_ordered.ff{data,index}
mv $outdir/hhsuite3_db_hhm_ordered.ffindex $outdir/hhsuite3_db_hhm.ffindex
mv $outdir/hhsuite3_db_hhm_ordered.ffdata $outdir/hhsuite3_db_hhm.ffdata
    
$hhblits_dir/ffindex_order sorting.dat $outdir/hhsuite3_db_a3m.ff{data,index} $outdir/hhsuite3_db_a3m_ordered.ff{data,index}
mv $outdir/hhsuite3_db_a3m_ordered.ffindex $outdir/hhsuite3_db_a3m.ffindex
mv $outdir/hhsuite3_db_a3m_ordered.ffdata $outdir/hhsuite3_db_a3m.ffdata



