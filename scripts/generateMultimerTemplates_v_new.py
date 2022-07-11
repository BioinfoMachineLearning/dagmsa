#this script takes input the fasta files of the templates according to the chains
#usage: python generateMultimerTemplates.py <input_fasta_directory> <output_dir> <--hhblits>
#1. Generates monomer templates-
#      (a) By aligning using hhblits, bfd and a custom monomer hhblits searchable database.
#      (b) Using the PIR alignment database.
#
#2. Generates an alignment file of all the templates in the PIR format. The alignments are separated by '/'
#3. Runs modeller on the file in (2)

import os, sys
from re import sub
from glob import glob

length_dict={} # length dictionary of each sequence in the alignment file especially the first one
chain_dict={}
s='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
chain2num_dict={}
reverse_chain2num_dict={}
#global MAX_ALN_NUM is the total number of subgraphs that will be generated
MAX_ALN_NUM=0
for i in range(len(s)):
    chain2num_dict[s[i]]=i+1
    reverse_chain2num_dict[i+1]=s[i]

fileNumber=0
print (chain2num_dict)
print (reverse_chain2num_dict)
#sys.exit()

def prob(e):
    return e[2]

# Define what this does: 


def readAlnPIRFile(alignment_pir_file):
    index_list=[]
    chain_alignment_dict={}
    first_line=True
    L=0
    first_fasta_index=""
    first_fasta=""
    with open (alignment_pir_file) as f:
        for line in f:
            if line.strip()=="": continue
            if line.startswith(">"):# and first_line: #The first line is a query or target
                query_index=line.strip().split(";")[1]
                line=f.readline()
                fasta=""
                while "*" not in line:
                    line=f.readline().strip()
                    fasta+=line
                #print ("L_F=",len(fasta)-1)
                chain_alignment_dict[query_index]=fasta
                #print (first_line)
                if first_line: 
                    if fasta.endswith("*"): 
                        L=len(fasta)-1
                    else:
                        L=len(fasta) 
                    first_fasta_index=query_index
                    first_fasta=fasta
                first_line=False
                index_list.append(query_index)
                fasta=""
            #elif line.startswith(">") and not first_line:
    return index_list, chain_alignment_dict, L, first_fasta, first_fasta_index

def removeAlns(infile,outfile,name):
    header="Ranked templates for "+name+"\n"
    contents=[]
    cont=[]
    i=0
    cont.append(header)
    with open (infile) as f:
        for line in f:
            if line.startswith("#"): continue
            #if line.startswith("1"): print("1111111 ONE 11111111111111")
            if ( line.strip().startswith("1") or line.strip().startswith("2") or line.strip().startswith("3") or line.strip().startswith("4") or line.strip().startswith("5") or line.strip().startswith("6") or line.strip().startswith("7") or line.strip().startswith("8") or line.strip().startswith("9") or line.strip().startswith("0")): 
                split=line.strip().split()
                #split[1]=split[1].upper()+"_"+str(i)
                i+=1
                new_line=""
                #contents.append([split[0],split[1],split[10],float(split[11])])
                contents.append([split[1],split[10],float(split[11])])
                #for i in range(len(split)-1):
                #    new_line+=split[i]+"\t"

                #contents.append(new_line.strip()+"\n")
    #print ("LENGTH_SPLIT #################=",len(split))
    #with open (infile) as f:
    #    c=f.readlines()
    #print ("LENGTH C = ",len(c))
    #print (cont)
    first=contents[0]
    contents.remove(first)
    #code for ranking the templates
    contents.sort(key=prob,reverse=False)
    contents.insert(0,first)
    k=1
    for item in contents:
        #cont.append(item[0]+"\t"+item[1]+"\t"+str(item[3])+"\t"+str(item[2])+"\n")
        cont.append(str(k)+"\t"+item[0]+"\t"+str(item[2])+"\t"+str(item[1])+"\n")
        k+=1
    with open (outfile,"w") as ff:
        ff.writelines(cont)


#read the index_list according to the chain number and select the most probable chains
def selectNonRedundantIndicesAccording2Chains(chain_index_list,chain_num,previous_chain_index_list):
    new_list=[]
    l=len(chain_index_list)
    i=0
    mini_common_list=[]
    while i<l-1:
        if chain_index_list[i][0:4]==chain_index_list[i+1][0:4]:
            print (chain_index_list[i],",",chain_index_list[i+1])
            mini_common_list.append(chain_index_list[i])
            #mini_common_list.append(chain_index_list[i+1])
            print ("i=",i,mini_common_list)
            if i==l-2 and len(mini_common_list)>=chain_num:
                mini_common_list.append(chain_index_list[i+1])
                print ("MCL=",mini_common_list,chain_num)
                new_list.append(mini_common_list[chain_num])
            i+=1
        else:
            if len(mini_common_list)>chain_num:
                print ("$$$$",len(mini_common_list),mini_common_list,chain_num,chain_index_list[i],chain_index_list[i+1])
                new_list.append(mini_common_list[chain_num])
            else:
                if previous_chain_index_list!=None:
                    if chain_index_list[i] not in previous_chain_index_list:
                        new_list.append(chain_index_list[i])
                else:
                    new_list.append(chain_index_list[i])

            mini_common_list=[]
            i+=1
            #pass
        if i==l-1:
            print ("MMMMM",i)
    
    #Last one gets missed.

    return new_list

def findUniqueNodes(all_indices):
    new_list=[]
    dag_dict={}
    for i in range(len(all_indices)):
        for chain_nodes in all_indices[i]:
            new_list.append(chain_nodes[0:4])
            if dag_dict.get(chain_nodes[0:4])==None:
                dag_dict[chain_nodes[0:4]]=[]
            dag_dict[chain_nodes[0:4]].append(chain_nodes)
    
    
    return sorted(list(set(new_list))), dag_dict

def addHypotheticalChains(nr_index_list):
    print (nr_index_list)
    chain_wise_node_list_nr=[]
    for nr_lst in nr_index_list:
        pass
    l=len(nr_index_list)
    i=0
    while i < l-1:
        first_node_list=nr_index_list[i]
        second_node_list=nr_index_list[i+1]



    return

def convertMSA2PIR(MSA_dict_, _first_index_list, _query_fasta_dict):
    file_contents=[]
    print ("Preparing PIR contents ...")
    MSAList=list(MSA_dict.keys())
    for i in range(len(MSAList)):
        row=MSA_dict_[MSAList[i]]
        #print (row)
        #print(type(row))
        master_key=MSAList[i]
        sub_list=row #list(row[master_key].keys())
        #print (sub_list)
        #print (type(sub_list))
        
        #sys.exit()
        #(">P1;"+rname[0:4]+"\n")
        #print (master_key, sub_list, type(sub_list))
        print (">P1;"+master_key+"\n")
        file_contents.append(">P1;"+master_key+"\n")
        tag="structure:"+master_key+"_"
        sequence=""
        for __i in range(len(sub_list)):
            chain_pdbs = list(sub_list[__i].keys())[0]
            tag+=chain_pdbs[4]
            sequence+=row[__i][chain_pdbs] #row[master_key][chain_pdbs]
            if sequence.endswith("*"): 
                sequence=sequence[:-1]
                sequence+="/"
        if sequence.endswith("/"):
            sequence=sequence[:-1]
            sequence+="*\n\n"
        tag+=": : : : : : : : \n"

        print (tag)
        print (sequence)
        file_contents.append(tag)
        file_contents.append(sequence)
        #print ("i=",i)

        #break
    
    #row=MSAList[-1]
    master_key=_first_index_list[0][0:4]
    sub_list=_first_index_list
    #print (_first_index_list)
    print(_query_fasta_dict)
    #sys.exit()
    print (">P1;"+master_key+"\n")
    file_contents.append(">P1;"+master_key+"\n")
    tag="sequence:"+master_key+"_"
    sequence=""
    #print ("Chain_pdbs=",chain_pdbs)
    #sys.exit()
    for chain_pdbs in sub_list:
        tag+=chain_pdbs[4]
        sequence+=_query_fasta_dict[chain_pdbs] #row[master_key][chain_pdbs]
        #print (sequence,"$$$$",chain_pdbs)
        if sequence.endswith("*"): 
            sequence=sequence[:-1]
            sequence+="/"
    #print (sequence)
    if sequence.endswith("/"):
        sequence=sequence[:-1]
        sequence+="*\n\n"
    tag+=": : : : : : : : \n"

    print (tag)
    print (sequence)
    file_contents.append(tag)
    file_contents.append(sequence)
    
    return file_contents


#########################################################################################3
# Main
##########################################################################################
inpdir=os.path.abspath(sys.argv[1])+"/" #Replace with a file which contains hard paths of all fastas
outdir=os.path.abspath(sys.argv[2])+"/"
hhblits_flag=False
if len(sys.argv)==4 and sys.argv[3]=="--hhblits":
    hhblits_flag=True
else:
    hhblits_flag=False

if not os.path.isdir(outdir): os.makedirs(outdir)
if not os.path.isdir(outdir+"fasta_pir/"): os.makedirs(outdir+"fasta_pir/")
if not os.path.isdir(outdir+"tmpdir/"): os.makedirs(outdir+"tmpdir/")
fastapirdir=outdir+"fasta_pir/"
tmpdir=outdir+"tmpdir/"

inp_pir_list=[] #the list of the input fasta files in pir format

chain_dict_list=[]

fastalist=sorted(glob(inpdir+"*.fasta"))
print (fastalist)
fileNumber=len(fastalist)
for file in fastalist:
    pir_name=fastapirdir+os.path.basename(file).replace(".fasta",".pir")
    print (pir_name)
    inp_pir_list.append(pir_name)
    #os.system("python fasta2pir.py "+file+" "+pir_name)

#Generate Alignment in PIR format using respective databases.
if not hhblits_flag:
    for file in inp_pir_list:
        name=os.path.basename(file).replace(".pir","")
        #os.system("python generateTemplateListPIR.py "+file+" "+tmpdir+os.path.basename(file).replace(".pir",""))
        removeAlns(tmpdir+name+"/"+name+"_build_profile.prf",tmpdir+name+"/"+name+".hhr",name)
        
else:
    pass
    """
    # Run the hhblits script for generating monomer a3m and ranking protocol
    for file in inp_pir_list:
        os.system("python generateTemplateListPIR.py "+file+" "+tmpdir+os.path.basename(file).replace(".pir",""))
    """

# Generate paired pir files,

d=[] #list of all the individual alignments in each alignment file as a dictionary
#Use d[i].keys() to get the key at the ith element
#d[0] is always the query

# Create the length dictionary of each fasta sequence from all the PIR alignment files. 
aln_file_dict={} #this stores all the alignment indices of each .pir alignment file

chain_wise_index_list=[]
chain_wise_length_list=[]
chain_wise_all_alignment_list=[]
first_index_list=[]
query_fasta_dict={}
chain_number=0

for file in inp_pir_list:
    name=os.path.basename(file).replace(".pir","")
    #print (tmpdir+name+"/"+name+"_aln.pir","77777777777777777")
    #dd,LL=createDict(tmpdir+name+"/"+name+"_aln.pir")
    index_list, chain_alignment_dict, L, first_fasta, first_fasta_index=readAlnPIRFile(tmpdir+name+"/"+name+"_aln.pir")
    index_list.remove(first_fasta_index)
    chain_wise_index_list.append(sorted(index_list))
    chain_wise_length_list.append(L)
    chain_wise_all_alignment_list.append(chain_alignment_dict)
    first_index_list.append(first_fasta_index)
    query_fasta_dict[first_fasta_index]=first_fasta

    #print(dd,"88888888888888888")
    #break
    #d.append(dd)
    #fn,L_nm=getLengthFromLengthDict(LL)
    #length_dict[fn]=L_nm

############################################################################################33

print (len(chain_wise_index_list))
chain_number=len(chain_wise_index_list)
#sys.exit()
#print (chain_wise_index_list)

#print (chain_wise_index_list[0])

nr_chain_wise_index_list=[]

for i in range(len(chain_wise_index_list)):
    if i == 0: 
        prev=None
    else:
        prev=[]
        for nr in nr_chain_wise_index_list:
            prev+=nr

    new_list_nr=selectNonRedundantIndicesAccording2Chains(chain_wise_index_list[i],i,prev)
    nr_chain_wise_index_list.append(new_list_nr)


#new_list_a=selectNonRedundantIndicesAccording2Chains(chain_wise_index_list[0],0,None)
#print ("NEW_A:",new_list_a, len(chain_wise_index_list[0]),len(new_list_a))
#print (sorted(chain_wise_index_list[0]))
#new_list_b=selectNonRedundantIndicesAccording2Chains(chain_wise_index_list[1],1,sorted(new_list_a))
#print ("NEW_B:",new_list_b, len(chain_wise_index_list[1]),len(new_list_b))

#new_list_c=selectNonRedundantIndicesAccording2Chains(chain_wise_index_list[2],2,new_list_a+new_list_b)
#print ("NEW_C:",new_list_c, len(chain_wise_index_list[2]),len(new_list_c))


#print (len(nr_chain_wise_index_list[0]))
#print (len(nr_chain_wise_index_list[1]))
#print (len(nr_chain_wise_index_list[2]))
#print (nr_chain_wise_index_list[1])

unique_nodes,dag_graph_index_list=findUniqueNodes(nr_chain_wise_index_list)

#hypo_nr_chain_wise_index_list=addHypotheticalChains(nr_chain_wise_index_list)

#sys.exit()

print (len(unique_nodes))
#unique_nodes.append(first_fasta_index[0:4])
print (unique_nodes)
print ("DAG=",dag_graph_index_list)
print("^^^^^^^^^^^^^^^^")
print (nr_chain_wise_index_list)

sorted_dag_keys=sorted(list(dag_graph_index_list.keys()))


print (sorted_dag_keys)


MSA_dict={}
for i in range(len(unique_nodes)):
    MSA_dict[unique_nodes[i]]=[]
    for chns in range(chain_number):
        MSA_dict[unique_nodes[i]].append({})

#print (MSA_dict)

for i in range(len(nr_chain_wise_index_list)):
    chain_wise_sub_list=nr_chain_wise_index_list[i]
    #print(chain_wise_sub_list)
    for chained_pdb in chain_wise_sub_list:
        pdb_code=chained_pdb[0:4]
        MSA_dict[pdb_code][i][chained_pdb]=chain_wise_all_alignment_list[i][chained_pdb]


    #print(chain_wise_sub_list)
for i in range(len(unique_nodes)):
    row=MSA_dict[unique_nodes[i]] #list of dictionaries
    #print (type(row))
    #print(row)
    #break
    pdb_name=unique_nodes[i]
    filled_chain_list=[]
    blank_chain_list=[]
    value_dict={}
    for chns in range(len(row)):
        if row[chns]=={}:
            blank_chain_list.append(chns)
            #if chn-1==
        else:
            filled_chain_list.append(chns)
            #print (row[chns])
            value_dict[chns]=list(row[chns].keys())[0][4]
    for i_blank in blank_chain_list:
        for j_filled in filled_chain_list:
            if i_blank<j_filled:
                j_chain=value_dict[j_filled]
                diff=j_filled-i_blank
                if chain2num_dict[j_chain]-diff not in filled_chain_list:
                    mod_chain=reverse_chain2num_dict[chain2num_dict[j_chain]-diff]
            if i_blank>j_filled:
                j_chain=value_dict[j_filled]
                diff=i_blank-j_filled
                if chain2num_dict[j_chain]+diff not in filled_chain_list:
                    mod_chain=reverse_chain2num_dict[chain2num_dict[j_chain]+diff]
        gaps="-"*chain_wise_length_list[i_blank]
        #if gaps
        MSA_dict[unique_nodes[i]][i_blank]={pdb_name+mod_chain:gaps+"*"}




print (MSA_dict)

#sys.exit()
#print (type(MSA_dict))
print (MSA_dict.keys())
file_contents=convertMSA2PIR(MSA_dict,first_index_list,query_fasta_dict)

if not os.path.isdir(outdir+"pir_out/"): os.makedirs(outdir+"pir_out/")
with open (outdir+"pir_out/"+name[0:4]+"_0.pir","w") as f:
    f.writelines(file_contents)
sys.exit()




#print (outdir+"pir_out/"+name)
#sys.exit()
if not os.path.isdir(outdir+"pir_out/"): os.makedirs(outdir+"pir_out/")
with open (outdir+"pir_out/"+name[0:4]+"_0.pir","w") as f:
    f.writelines(file_contents)

sys.exit()
dag_aln_list=[]
for dag_k in sorted_dag_keys:
    sub_list=sorted(dag_graph_index_list[dag_k])
    print (sub_list)
    d_aln={}
    for _i in range(len(sub_list)):
        sub_list_chain=sub_list[_i]
        if chain_wise_all_alignment_list[_i].get(sub_list_chain)!=None:
            d_aln[sub_list_chain]=chain_wise_all_alignment_list[_i][sub_list_chain]
        else:
            current_chain=sub_list_chain[4]
            chain_code=reverse_chain2num_dict[chain2num_dict[current_chain]-1]
            d_aln[sub_list_chain[0:4]+chain_code]="-"*chain_wise_length_list[_i]
    dag_aln_list.append(d_aln)

#print (dag_aln_list[9])
#print (sorted_dag_keys)
        



