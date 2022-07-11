from doctest import master
from itertools import chain
import os, sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from numpy import rank

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
        if chain_num < len(sublist):             
            print ("chain_num=",chain_num,sublist)
            selectedKeys.append(sublist[chain_num])
    
    return sorted(selectedKeys), pdb_code_keys

def prob(e):
    return e[-1]


def rank2dict(rank_list):
    l=[]
    l_d={}
    #print (type(rank_list))
    #print(rank_list[0])

    #print (len(rank_list))
    for line in rank_list:
        split=line.strip().split()
        l.append([int(split[0]),split[1],float(split[2]),float(split[3])])
        l_d[split[1]]=[]
        l_d[split[1]].append(float(split[2]))
        l_d[split[1]].append(float(split[3]))
    
    #print (l)

    return l,l_d

def generateOverAllRanks(selected_k_dict,rfl,all_ranks):
    
    master_rank=[]
    print (type(selected_k_dict))
    print (selected_k_dict)
    for k, items in selected_k_dict.items():
        print (k, len(items))

    print (rfl)
    print (len(all_ranks))
    reformated_rank=[]
    for i in range(len(all_ranks)):
        rl,rank_d=rank2dict(all_ranks[i])
        reformated_rank.append(rank_d)
    
    selected_pdb_keys=list(selected_k_dict.keys())
    for k in selected_pdb_keys:
        chainwise_key_sublist=selected_k_dict[k]
        eval=0
        seqsim=0
        target=""
        chns=""
        for _k_ in range(len(chainwise_key_sublist)):
            print ("_k_=",_k_,chainwise_key_sublist[_k_], len(chainwise_key_sublist[_k_]))
            chain_key=chainwise_key_sublist[_k_]
            if len(chainwise_key_sublist) == len(reformated_rank):
                if reformated_rank[_k_].get(chain_key)!=None: 
                    target=chain_key
                    eval+=reformated_rank[_k_][chain_key][0]
                    seqsim+=reformated_rank[_k_][chain_key][1]
                print (":",target,eval,seqsim,":")
                chns+=target[4]
                #print ("chain_key=",chain_key, "reformated_rank=",reformated_rank[_k_][chain_key])
                #get ranks from chainwise_key[x] and reformated_rank[x]
                pass
            else: #this is the issue with 
                #while 
                print ("ELSE:","_k_=",_k_,chainwise_key_sublist[_k_], len(chainwise_key_sublist[_k_]))
                if reformated_rank[_k_].get(chain_key)==None:
                    target=chain_key
                    chns+=reverse_chain2num_dict[chain2num_dict[target[4]]+1]
                    eval+=reformated_rank[_k_][chain_key][0]
                    seqsim+=reformated_rank[_k_][chain_key][1]
                    print ("HERE:::",chns,"target=",target)
                    #chns 
                    #print ("OTHER_CHAIN_KEYS=",chain_key)
                    continue
                else:
                    print ("OTHER_CHAIN_KEYS=",chain_key)
                    target=chain_key
                    chns+=target[4]
                    eval+=reformated_rank[_k_][chain_key][0]
                    seqsim+=reformated_rank[_k_][chain_key][1]
                    print ("HERE:::",chns,"target=",target)
                #ref_rank=reformated_rank[_k_][chain_key]
                #print ("chain_key=",chain_key, "reformated_rank=",reformated_rank[_k_][chain_key])
                pass
        master_rank.append([target[0:4]+"_"+chns,eval,seqsim])
        #for chain_key in chainwise_key_sublist:
            #for d in rank_d:

    print (master_rank)


    return

monomer_pir_file_list=[]
monomer_pir_file_list_file=os.path.abspath(sys.argv[1])
chain_number=int(sys.argv[2])

selected_keys_dict={}

with open (monomer_pir_file_list_file) as f:
    for line in f:
        monomer_pir_file_list.append(line.strip())

all_sorted_rank_list=[]
rank_file_list=[]
for ii in range(len(monomer_pir_file_list)):
    monomer_pir_file=monomer_pir_file_list[ii] #os.path.abspath(sys.argv[1])
    print ("working on ",monomer_pir_file)
    rank_file=os.path.dirname(monomer_pir_file)+"/"+os.path.basename(monomer_pir_file).replace("_aln.pir",".hhr")
    print ("RANK_FILE=",rank_file)
    

    outdir=os.path.dirname(monomer_pir_file)+"/patch_monomer/"
    if not os.path.isdir(outdir): os.makedirs(outdir)

    chain_name_list,align_list=readPIRFile(monomer_pir_file)
    first_alignment_dict={chain_name_list[0]:align_list[0]}
    chained_alingment_dict={} #The first sequence is the query. This is ignored here. 
    blanks="-"*(len(align_list[0]))
    for j in range(1,len(chain_name_list)):
        chained_alingment_dict[chain_name_list[j]]=align_list[j]


    #print
    # The selectedKeys are the NR keys, pdb_code_list is the 4 leter PDB code 
    selectedKeys, pdb_code_list=getNRAlignmentIndex(chained_alingment_dict, ii, blanks) #getNRAlignmentIndex(chained_alingment_dict, chain_number-1, blanks)

    print ("Selected_KEYS=",selectedKeys)

    for skeys in selectedKeys:
        if selected_keys_dict.get(skeys[0:4])!=None:
            selected_keys_dict[skeys[0:4]].append(skeys)
        else:
            selected_keys_dict[skeys[0:4]]=[]
            selected_keys_dict[skeys[0:4]].append(skeys)


    

    rank_dict_contents={}
    #i=1
    lab=""
    with open (rank_file) as f:
        for line in f:
            if line.strip().startswith("Ranked"): 
                lab=line
                continue
            split=line.strip().split()
            nm=line.strip().split()[1].split("_")[0]
            #print ("split[1].split("")[0]=",split[1].split("_")[0])
            rank_dict_contents[nm] = split[1].split("_")[0]+"\t"+split[2]+"\t"+split[3]+"\n"
            #print ("rank_dict_contents[nm]=",rank_dict_contents[nm])

    nr_rank=[]
    #nr_rank.append("Ranked templates for "+list(first_alignment_dict.keys())[0]+"\n")
    i=0
    for nr_keys in selectedKeys:
        i+=1
        #split=
        #print ("i==",str(i))
        nr_rank.append([rank_dict_contents[nr_keys],float(rank_dict_contents[nr_keys].split()[-1])])
        #print ("BR_RANK=",nr_rank[-1])
        #print(rank_dict_contents[nr_keys].split()[-1])

    #print (nr_rank[0])
    nr_rank.sort(key=prob,reverse=True)
    sorted_nr_rank=nr_rank
    #print (type(sorted_nr_rank))
    sorted_rank_list=[]
    i_count=0
    for itm in sorted_nr_rank:
        i_count+=1
        #print ("I_count=",i_count,"itm=",itm)
        sorted_rank_list.append(str(i_count)+"\t"+itm[0])

    #break
    all_sorted_rank_list.append(sorted_rank_list.copy())
    sorted_rank_list.insert(0,lab)

    with open (outdir+os.path.basename(monomer_pir_file).replace("_aln.pir",".hhr"),"w") as f:
        f.writelines(sorted_rank_list)

    rank_file_list.append(outdir+os.path.basename(monomer_pir_file).replace("_aln.pir",".hhr"))

    if (ii==len(monomer_pir_file_list)-1): print ("Selected_Keys_Dict=",selected_keys_dict)   

overall_rank=generateOverAllRanks(selected_keys_dict,rank_file_list,all_sorted_rank_list)


    #print(align_list)
    #print (len(align_list))
    #print ()

#print ("THIS=",all_sorted_rank_list[2])
#print (len(all_sorted_rank_list))
sys.exit()

for _i in range(len(all_sorted_rank_list)-1):
    for _j in range(_i,len(all_sorted_rank_list)):
        split_A=all_sorted_rank_list[_i].split()
        split_B=all_sorted_rank_list[_j].split()



final_chain_keys={}

for k,value in selected_keys_dict.items():
    for _i in range(len(value)):

        pass
    #if len(value)!=


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
        #print ("Orig Fastas A")
        #print (fasta_A)
        #print ("Template Fastas A")
        #print (fasta_template_A)
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
        #print ("########################3")

        #print (best_aln_B[0])
        #print (best_aln_B[1])
        #orig_aln.append(best_aln_A[0]+"/"+best_aln_B[0]+"*")
        #local_aln.append(best_aln_A[1]+"/"+best_aln_B[1]+"*")
        after_realign.append([best_aln_A[0],best_aln_A[1]]) #first is the query

        #print (len(best_aln_A[0]),len(best_aln_A[1]))
        #sys.exit()
        """
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
        """

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
