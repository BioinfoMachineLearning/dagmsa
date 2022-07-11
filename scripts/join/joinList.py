import os, sys


a3mfolder="./single_a3ms/"
outfolder="./joined_a3ms/"
if not os.path.isdir(outfolder): os.makedirs(outfolder)

l=[]

with open ("lll.txt") as f:
    for line in f:
        l.append(line.strip())
        
        
for pair in l:
    A=pair.split("_")[0]
    B=pair.split("_")[1]
    os.system("python joinAlignments4Hetero_v2.py "+a3mfolder+A+".a3m "+a3mfolder+B+".a3m "+outfolder)
    