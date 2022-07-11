from itertools import chain
import os, sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

s='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
s+=s+s.lower()
chain2num_dict={}
reverse_chain2num_dict={}
#global MAX_ALN_NUM is the total number of subgraphs that will be generated
MAX_ALN_NUM=0
for i in range(len(s)):
    chain2num_dict[s[i]]=i+1
    reverse_chain2num_dict[i+1]=s[i]


def readFastaFile(file):
    fasta=""
    with open (file) as f:
        for line in f:
            if line.startswith(">"): continue
            fasta+=line.strip()
    return fasta

def getMaxConsecutiveGaps(sequence):
    num_of_gap_sequence=0
    if ("-" in sequence):
        gap_found=False
        i=-1
        while i < len(sequence)-1:
            i+=1
            if (sequence[i]=="-"):
                gap_found=True
                for j in range(i,len(sequence)):
                    
                    if (sequence[j]!="-"):
                        num_of_gap_sequence+=1
                        i=j-1
                        break
        #pass
    else:
        return 0
            
    return num_of_gap_sequence   

def selectBestAlignment(alignment):
    best=getMaxConsecutiveGaps(alignment[0][0])+getMaxConsecutiveGaps(alignment[0][1])
    best_aln=alignment[0]
    for i in range(1,len(alignment)):
        score_a1=getMaxConsecutiveGaps(alignment[i][0])
        score_a2=getMaxConsecutiveGaps(alignment[i][1])
        #print (score_a1,score_a2)
        if (score_a1+score_a2<=best):
            best=score_a1+score_a2
            best_aln=alignment[i]
        
    return best_aln


def readLocalFile(file,name,orig_fasta,rank_list):
    local_aln={}
    orig_aln={}
    with open (file) as f:
        for line in f:
            if line.strip()=="": continue
            for rname in rank_list:
                if rname in line:
                    f.readline()
                    line=f.readline()
                    nogap_seq=line.strip()#.replace("-","")
                    alignments = pairwise2.align.globalms(orig_fasta, nogap_seq,5,-4,-1,-0.1)
                    best_aln=selectBestAlignment(alignments)
                    local_aln[rname]=best_aln[1]+"*"
                    orig_aln[rname]=best_aln[0]+"*"
                    #print (len(best_aln[1].strip()))

        

    return local_aln, orig_aln

def readPIRFile(filename):
    aln_list=[]
    name_list=[]
    with open (filename) as f:
        for line in f:
            if line.startswith(">"):
                name=line.strip().split(";")[1]
                name_list.append(name)
                line=f.readline()
                fasta=""
                while "*" not in line:
                    line=f.readline().strip()
                    fasta+=line
                #print ("L_F=",len(fasta)-1)
                fasta=fasta.replace("*","")
                aln_list.append(fasta)
                #print (first_line)
                
    return name_list, aln_list


def getNRAlignmentIndex(chained_alingment_dict, chain_num, bln):
    all_keys=sorted(list(chained_alingment_dict.keys()))
    pdb_code_keys=[]
    for key in all_keys:
        pdb_code_keys.append(key[0:4])
    pdb_code_keys=sorted(list(set(pdb_code_keys)))
    print ("PDBS=",pdb_code_keys)
    selectedKeys=[]
    key_found=False
    for pdb_code in pdb_code_keys:
        sublist=[]
        for key in all_keys:
            if pdb_code==key[0:4]:
                sublist.append(key)
                key_found=True
            else:
                if key_found==True: 
                    key_found=False
                    break
        if chain_num <= len(sublist):             
            print (sublist,chain_num)
            selectedKeys.append(sublist[chain_num])
    
    return sorted(selectedKeys), pdb_code_keys

def prob(e):
    return e[-1]




monomer_pir_file=os.path.abspath(sys.argv[1])
rank_file=os.path.dirname(monomer_pir_file)+"/"+os.path.basename(monomer_pir_file).replace("_aln.pir",".hhr")
chain_number=int(sys.argv[2])

outdir=os.path.dirname(monomer_pir_file)+"/patch_monomer/"
if not os.path.isdir(outdir): os.makedirs(outdir)

chain_name_list,align_list=readPIRFile(monomer_pir_file)
first_alignment_dict={chain_name_list[0]:align_list[0]}
chained_alingment_dict={} #The first sequence is the query. This is ignored here. 
blanks="-"*(len(align_list[0]))
for i in range(1,len(chain_name_list)):
    chained_alingment_dict[chain_name_list[i]]=align_list[i]


print 
selectedKeys, pdb_code_list=getNRAlignmentIndex(chained_alingment_dict, chain_number-1, blanks)

print (selectedKeys)


rank_dict_contents={}
i=1
lab=""
with open (rank_file) as f:
    for line in f:
        if line.strip().startswith("Ranked"): 
            lab=line
            continue
        split=line.strip().split()
        nm=line.strip().split()[1].split("_")[0]
        rank_dict_contents[nm] = split[1].split("_")[0]+"\t"+split[2]+"\t"+split[3]+"\n"

nr_rank=[]
#nr_rank.append("Ranked templates for "+list(first_alignment_dict.keys())[0]+"\n")
i=0
for nr_keys in selectedKeys:
    i+=1
    #split=
    nr_rank.append([str(i)+"\t"+rank_dict_contents[nr_keys],float(rank_dict_contents[nr_keys].split()[-1])])
    #print(rank_dict_contents[nr_keys].split()[-1])

#print (nr_rank[0])
nr_rank.sort(key=prob,reverse=True)
sorted_nr_rank=nr_rank
#print (type(sorted_nr_rank))
sorted_rank_list=[]
for itm in sorted_nr_rank:
    sorted_rank_list.append(itm[0])

sorted_rank_list.insert(0,lab)

with open (outdir+os.path.basename(monomer_pir_file).replace("_aln.pir",".hhr"),"w") as f:
    f.writelines(sorted_rank_list)


        



#print(align_list)
#print (len(align_list))
#print ()


#sys.exit()





#fasta_file_A=os.path.abspath("/data/farhan/multicom4s/tool/multicom4s_TMB/scripts/homotrimer/homotrimer_fasta/1RERA.fasta")
#fasta_file_B=os.path.abspath("/data/farhan/multicom4s/tool/multicom4s_TMB/scripts/homotrimer/homotrimer_fasta/1RERB.fasta")
#output_folder=os.path.abspath("./out")+"/"
#template_fasta_folder=os.path.abspath("./")
#if not os.path.isdir(output_folder): os.makedirs(output_folder)

local_aln=[]
orig_aln=[]

after_realign=[]

fasta_A= align_list[0] #readFastaFile(fasta_file_A)

for i in range(len(selectedKeys)):
    #fasta_B=align_list[i] #readFastaFile(fasta_file_B)
    #print (os.path.basename(fasta_file_A))
    name=selectedKeys[i] #os.path.basename(fasta_file_A).replace(".fasta","")+"_"+os.path.basename(fasta_file_B).replace(".fasta","")

    rname=name #"1LD4A_1LD4B"

    template_name_A=rname.split("_")[0]
    #template_name_B=rname.split("_")[1]
    fasta_template_A=chained_alingment_dict[selectedKeys[i]] #readFastaFile(template_fasta_folder+"/"+template_name_A+".fasta")
    #fasta_template_B=readFastaFile(template_fasta_folder+"/"+template_name_B+".fasta")
    print ("Orig Fastas A")
    print (fasta_A)
    print ("Template Fastas A")
    print (fasta_template_A)
    #print ("Orig Fastas B")
    #print (fasta_B)
    #print ("Template Fastas B")
    #print (fasta_template_B)
    alignments_A=pairwise2.align.globalms(fasta_A, fasta_template_A,5,-4,-1,-0.1)
    #alignments_B=pairwise2.align.globalms(fasta_B, fasta_template_B,5,-4,-1,-0.1)
    best_aln_A=selectBestAlignment(alignments_A)
    #best_aln_B=selectBestAlignment(alignments_B)

    #print (best_aln_A[0])
    #print (best_aln_A[1])
    print ("########################3")

    #print (best_aln_B[0])
    #print (best_aln_B[1])
    #orig_aln.append(best_aln_A[0]+"/"+best_aln_B[0]+"*")
    #local_aln.append(best_aln_A[1]+"/"+best_aln_B[1]+"*")
    after_realign.append([best_aln_A[0],best_aln_A[1]]) #first is the query

    print (len(best_aln_A[0]),len(best_aln_A[1]))
    #sys.exit()
    #"""
    with open (output_folder+"/"+name+"_"+str(0)+".pir","w+") as f:
            f.write(">P1;"+rname[0:5]+"\n")
            chain_A=rname[4]
            #chain_A=rname.strip().split("_")[0][-1]
            #chain_B=rname.strip().split("_")[1][-1]
            #f.write("structureN:"+rname[0:4]+"::"+chain_A+"::"+chain_B+"::::\n")
            f.write("structure:"+rname[0:5]+"::"+chain_A+"::"+""+"::::\n")
    #        f.write(local_aln[i]+"\n")
            f.write([0]+"\n")
            f.write("\n")
            f.write(">P1;"+name+"\n")
            targ_chain_A=name.strip().split("_")[0][-1]
            targ_chain_B=name.strip().split("_")[1][-1]
            f.write("sequence: :"+targ_chain_A+" : :"+targ_chain_B+" : : : : : \n")
    #        f.write(orig_aln[i]+"\n")
            f.write(orig_aln[0]+"\n")
    #"""

"""
with open (outdir+os.path.basename(monomer_pir_file),"w+") as f:
            f.write(">P1;"+rname[0:5]+"\n")
            chain_A=rname[4]
            #chain_A=rname.strip().split("_")[0][-1]
            #chain_B=rname.strip().split("_")[1][-1]
            #f.write("structureN:"+rname[0:4]+"::"+chain_A+"::"+chain_B+"::::\n")
            f.write("structure:"+rname[0:5]+"::"+chain_A+"::"+""+"::::\n")
    #        f.write(local_aln[i]+"\n")
            f.write([0]+"\n")
            f.write("\n")
            f.write(">P1;"+name+"\n")
            targ_chain_A=name.strip().split("_")[0][-1]
            targ_chain_B=name.strip().split("_")[1][-1]
            f.write("sequence: :"+targ_chain_A+" : :"+targ_chain_B+" : : : : : \n")
    #        f.write(orig_aln[i]+"\n")
            f.write(orig_aln[0]+"\n")
"""
