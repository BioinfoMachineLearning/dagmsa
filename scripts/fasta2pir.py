#this script takes input a fasta file, converts it to .pir file, 
#usage: python fasta2pir.py <inputfasta> <outputpirfile>
from modeller import *
import os, sys

fasta=os.path.abspath(sys.argv[1])
outfile=os.path.abspath(sys.argv[2])
name=os.path.basename(fasta).replace(".fasta","")

env = Environ()
aln_inp = Alignment(env)
aln_inp.append(file=fasta,alignment_format="FASTA",align_codes="ALL")
aln_inp.write(file=outfile, alignment_format='PIR')