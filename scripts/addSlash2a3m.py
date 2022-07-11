# this script adds a "/" after the first fasta sequence to separate the two 
# fasta sequences in the alignment a3m file

import os, sys
a3mfile=os.path.abspath(sys.argv[1])
fasta_file_A=os.path.abspath(sys.argv[2])
fasta_file_B=os.path.abspath(sys.argv[3])
outputfile=os.path.abspath(sys.argv[4])

def readFastaFile(file):
    fasta=""
    with open (file) as f:
        for line in f:
            if line.startswith(">"): continue
            fasta+=line.strip()
    return fasta

len_A=len(readFastaFile(fasta_file_A))
len_B=len(readFastaFile(fasta_file_B))

a3m=[]
with open (a3mfile) as f:
    for line in f:
        if line.startswith(">"):
            a3m.append(line)
        else:
            line=line[0:len_A]+"/"+line[len_A:]
            a3m.append(line)

with open (outputfile,"w") as f:
    f.writelines(a3m)


