#this script takes input a fasta file, converts it to .pir file, and then generates a list of the templates
#usage: python generateTemplateList.py <inputfasta> <output folder>

from modeller import *
#from Bio import SeqIO
import os, sys

fasta=os.path.abspath(sys.argv[1])
outfolder=os.path.abspath(sys.argv[2])+"/"
name=os.path.basename(fasta).replace(".fasta","")
if not os.path.isdir(outfolder): os.makedirs(outfolder)


#records = SeqIO.parse(fasta, "fasta")
#count = SeqIO.write(records, outfolder+name+"_inp.pir", "pir")
#print("Converted %i records" % count)

#sys.exit()

if fasta.endswith(".fasta"):
    env = Environ()
    aln_inp = Alignment(env)
    aln_inp.append(file=fasta,alignment_format="FASTA",align_codes="ALL")
    aln_inp.write(file=outfolder+name+"_inp.pir", alignment_format='PIR')
    name+="_inp.pir"


print("########################NAME=",name,"  OUTDIR=",outfolder)

log.verbose()
env = Environ()
#sys.exit()
# Read in the sequence database in binary format
sdb = SequenceDB(env, seq_database_file='/data/farhan/multicom4s/tool/multicom4s_TMB/databases/multimer_pir_db.hdf5',#'/data/farhan/multicom4s/tool/multicom4s_TMB/databases/pir/pdball.hdf5',
                 seq_database_format='BINARY', chains_list='ALL')

#sys.exit()

# Read in the target sequence in PIR alignment format
aln = Alignment(env)
#print (aln)
#sys.exit()
if not fasta.endswith(".fasta"):
    print("$$$$$$$$$$$$$$$$$$$$$  "+outfolder+name)
    aln.append(file=fasta, alignment_format='PIR', align_codes='ALL')
else:
    print("$$$$$$$$$$$$$$$$$$$$$  "+outfolder+name)
    aln.append(file=outfolder+name, alignment_format='PIR', align_codes='ALL')

# Convert the input sequence "alignment" into profile format
prf = aln.to_profile()

# Scan sequence database to pick up homologous sequences
prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
          gap_penalties_1d=(-500, -50), n_prof_iterations=1,
          check_profile=False, max_aln_evalue=0.01)

# Write out the profile in text format
prf.write(file=outfolder+name.replace(".pir","_")+'build_profile.prf', profile_format='TEXT')

# Convert the profile back to alignment format
aln = prf.to_alignment()

#print (aln)
#print (type(aln))
#print (len(aln))
#print (aln[0][0])
#sys.exit()

#- Write out a PIR alignment file
#if name.endswith(".pir"): name.replace(".pir","")
aln.write(file=outfolder+name.replace(".pir","")+'_aln.pir', alignment_format='PIR')
