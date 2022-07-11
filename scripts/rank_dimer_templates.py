#this script will combine the ranked monomer templates into dimer template list
#the ranking is recalculated using average

import os, sys

template_A="/data/farhan/multicom4s/tool/multicom4s_TMB/work_dir/5GNI_out/5GNIA.rank"
template_B="/data/farhan/multicom4s/tool/multicom4s_TMB/work_dir/5GNI_out/5GNIB.rank"

dimer_template="/data/farhan/multicom4s/tool/multicom4s_TMB/work_dir/5GNI_out/5GNI_AB.rank"

temp_A=[]
temp_B=[]
header=""
with open (template_A) as f:
    for line in f:
        if line.startswith("Ranked"): 
            header=line
            continue
        temp_A.append(line.strip())

with open (template_B) as f:
    for line in f:
        if line.startswith("Ranked"): continue
        temp_B.append(line.strip())
        
print (len(temp_A)," templates found for monomer ",os.path.basename(template_A).split(".")[0])
print (len(temp_B)," templates found for monomer ",os.path.basename(template_B).split(".")[0])

temp_AB=[]
#temp_AB.append("Rank\tName\tE-Value\tProb\n")

for A in temp_A:
    split=A.split()
    rank_A=split[0]
    name_A=split[1]
    eval_A=split[2]
    prob_A=split[3]
    #print (split)
    for B in temp_B:
        split_B=B.split()
        rank_B=split_B[0]
        name_B=split_B[1]
        eval_B=split_B[2]
        prob_B=split_B[3]
        if name_A==name_B: continue #same structure found as template
        if name_A[0:4]==name_B[0:4]:
            name_AB=name_A+","+name_B
            eval_AB=(float(eval_A)+float(eval_B))/2
            prob_AB=(float(prob_A)+float(prob_B))/2
            temp_AB.append([name_AB,eval_AB,prob_AB])
        
        
    #break
print (len(temp_AB))

with open ("dimer_list.txt","w") as f:
    for AB in temp_AB:
        f.write(AB[0]+"\n")







