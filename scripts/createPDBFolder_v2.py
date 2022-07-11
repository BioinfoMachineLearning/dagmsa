#this script searches the atom folder. copies the .atm files and merges them into one pdb file according to the rank_file

import os, sys

from glob import glob

#from numpy import rank


rankfile=sys.argv[1]
atom_dir=os.path.abspath(sys.argv[2])+"/"
outdir=os.path.abspath(sys.argv[3])+"/"

if not os.path.isdir(outdir):os.makedirs(outdir)

ranks=rankfile
i=0
"""
with open (rankfile) as f:
    for line in f:
        if line.startswith("Ranked"): continue
        ranks.append(line.strip().split()[1])
        i+=1
        #if i>=10: break
print (len(ranks))
"""

pair=glob(atom_dir+ranks+"*.atm")
print (ranks)


print (pair)
#sys.exit()
contents=[]
split=pair
pdb_file_name=ranks+".pdb"
for atom_file in sorted(split):
    #print (atom_file)
    with open (atom_file) as f:
        for line in f:
            #if line.startswith("ATOM"):
            #    if line[21]=="" or line[21]==" ":
                    #print (line)
            #        line=line[0:21]+os.path.basename(atom_file)[4]+line[22:]
                    #print ("New line=",line)
            if line.startswith("END\n"): 
                line=line.replace("END\n","TER\n")
                contents.append(line)
                continue
            if line.startswith("END"): line=line.replace("END","TER\n")
            contents.append(line)
contents.append("END\n")
with open (outdir+pdb_file_name,"w") as ff:
    ff.writelines(contents)

#print (contents)
#print (atom_file)
#print (split)