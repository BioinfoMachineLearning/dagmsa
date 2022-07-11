from turtle import pd


file="/data/farhan/multicom4s/tool/multicom4s_TMB/databases/pir/pdball.pir"

atom_list=[]
pdb_list=[]

with open(file) as f:
    for line in f:
        if line.startswith(">"):
            split=line.strip().split(";")[1]
            atom_list.append(split+"\n")
            pdb_list.append(split[0:4]+"\n")

pdb_list=list(set(pdb_list))

with open ("pdball_list.txt","w") as f:
    f.writelines(pdb_list)

with open ("pdball_atom_list.txt","w") as f:
    f.writelines(atom_list)
