#this script searches the atom folder. copies the .atm files and merges them into one pdb file according to the rank_file

import os, sys

from glob import glob


rankfile=os.path.abspath(sys.argv[1])
atom_dir=os.path.abspath(sys.argv[2])+"/"
outdir=os.path.abspath(sys.argv[3])+"/"

if not os.path.isdir(outdir):os.makedirs(outdir)

ranks=[]
i=0
with open (rankfile) as f:
    for line in f:
        if line.startswith("Ranked"): continue
        ranks.append(line.strip().split()[1])
        i+=1
        if i>=10: break
print (len(ranks))

for pair in ranks:
    print (pair)
    contents=[]
    pdb_file_name=pair[0:4]+".pdb"
    split=pair.split("_")
    for atom_file in split:
        with open (atom_dir+atom_file+".atm") as f:
            for line in f:
                if line.startswith("END\n"): 
                    line=line.replace("END\n","TER\n")
                    contents.append(line)
                    continue
                if line.startswith("END"): line=line.replace("END","TER\n")
                contents.append(line)
    contents.append("END\n")
    with open (outdir+pdb_file_name,"w") as ff:
        ff.writelines(contents)



