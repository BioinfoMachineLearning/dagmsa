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
                split[1]=split[1].upper()+"_"+str(i)
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

def createDict(aln_pir_file): 
    #d={}
    d=[]
    i=-1
    L_chains=[]
    fasta=""
    idx=""
    idx_prev=""
    global MAX_ALN_NUM
    with open (aln_pir_file) as f:
        for line in f:
            if line.strip()=="": continue
            if line.startswith(">"):
                #if len(d)!=0: d.append(dd)
                #dd=[]
                i+=1
                #idx=line.strip().split()[0].split(";")[1].upper()+"_"+str(i)
                idx=line.strip().split()[0].split(";")[1]
                
                line=f.readline() #sequence:...
                if fasta=="": # 1st sequence
                    idx_prev=idx
                    
                if fasta!="": # 1st sequence                    
                    d.append({idx_prev:fasta})                    
                    #length_dict
                    L_chains.append({idx_prev:len(fasta)})
                    fasta=""
                    #i+=1
                    idx_prev=idx
                
                #continue
            else:
                fasta+=line.strip()
        d.append({idx:fasta})
        L_chains.append({idx:len(fasta)})
    
    #print (len(d[0]))
    #print (d[0])
    print (len(d))
    
    m=0
    for k in d:
        #keys=k.keys()
        m=len(k)
        #print ("M=",m)
        if MAX_ALN_NUM<m: MAX_ALN_NUM=m
    """
    for k in d:
        keys=k.keys()
        m=len(keys)
        if MAX_LENGTH<m: MAX_LENGTH=m
    """
    return d, L_chains
                
#####
def createAlnDict(old_dict,index_list,aln_pir_file): 
    #d={}
    d=[]
    i=-1
    L_chains=[]
    fasta=""
    idx=""
    idx_prev=""
    global MAX_ALN_NUM
    with open (aln_pir_file) as f:
        for line in f:
            if line.strip()=="": continue
            if line.startswith(">"):
                #if len(d)!=0: d.append(dd)
                #dd=[]
                i+=1
                #idx=line.strip().split()[0].split(";")[1].upper()+"_"+str(i)
                idx=line.strip().split()[0].split(";")[1]
                
                line=f.readline() #sequence:...
                if fasta=="": # 1st sequence
                    idx_prev=idx
                    
                if fasta!="": # 1st sequence                    
                    d.append({idx_prev:fasta})                    
                    #length_dict
                    L_chains.append({idx_prev:len(fasta)})
                    fasta=""
                    #i+=1
                    idx_prev=idx
                
                #continue
            else:
                fasta+=line.strip()
        d.append({idx:fasta})
        L_chains.append({idx:len(fasta)})
    
    #print (len(d[0]))
    #print (d[0])
    print (len(d))
    
    m=0
    for k in d:
        #keys=k.keys()
        m=len(k)
        #print ("M=",m)
        if MAX_ALN_NUM<m: MAX_ALN_NUM=m
    """
    for k in d:
        keys=k.keys()
        m=len(keys)
        if MAX_LENGTH<m: MAX_LENGTH=m
    """
    return d, L_chains

#####




def findInterface(code, targetDict):
    d={}
    keys=list(targetDict.keys()) #sorted list ensures chains are in order
    print ("CHABIS=",len(keys))
    return d

def joinInterface(code, targetDict):
    return

def lengthOfAlignedFasta(d):
    for key,value in d.items():
        return len(value)
    return None

def createGraph(d):
    sub_graphs=[] #number of subgraphs will be = MAX_LENGTH
    for i in range(len(d)):
        k = d[i].keys()
        node_dict={}
        for key in k:
            node_dict[key[0:4]]=[]


    node_list=[]

    return

def generateChainDict(d): #This generates the 
    #dictionary of lists
    chain_list=[]
    #ch_dict={}
    i=0
    for aln_dict in d:
        ch_dict={}
        chain_list.append(ch_dict)
        keys = aln_dict.keys()
        for k in keys:
            k_pdb=k[0:4]
            if ch.get(k_pdb)==None:
                ch_dict[k_pdb]=[]
                ch_dict[k_pdb].append(k[4])
                #ch_dict[k_pdb].append(k)
#                pass
            else:
                ch_dict[k_pdb].append(k[4])
                #ch_dict[k_pdb].append(k)
                #pass

    return ch_dict

"""
def setUpChainwiseLengths(ld_list):
    L_dict={}
    for i in range(len(ld_list)):


    return
"""
####Problem HERE
"""
def removeRedundants(ld1,ld2): #### Conditions ####problem removing from ld1, ld2 but still existing in old
    new_ld1=[]
    new_ld2=[]
    old_ld1=ld1.copy()
    old_ld2=ld2.copy()
    for i in range(len(old_ld1)):
        for j in range(len(old_ld2)):
            if old_ld1[i][0:4]==old_ld2[j][0:4]:
                if old_ld1[i][4]==old_ld2[j][4]:
                    ld2.remove(old_ld2[j])
                    continue
                if (chain2num_dict[old_ld2[j][4]]-chain2num_dict[old_ld1[i][4]]==1):
                    if old_ld2[j] in ld1: ld1.remove(old_ld2[j])
                    continue
                #if (abs(chain2num_dict[old_ld2[j][4]]-chain2num_dict[old_ld1[i][4]])>1):
                #if (abs(chain2num_dict[old_ld1[i][4]]-chain2num_dict[old_ld2[j][4]])>1):
                #if (chain2num_dict[old_ld1[i][4]]-chain2num_dict[old_ld2[j][4]]<1):
                #if (old_ld1[i][4]==old_ld2[j][4]) or abs((chain2num_dict[old_ld1[i][4]]-chain2num_dict[old_ld2[j][4]])>1):
                    #print (old_ld1[i],",",old_ld2[j],"((((((((((((((((((")
                    #if old_ld1[i] in ld1: ld1.remove(old_ld1[i])
                    #if old_ld2[j] in ld2: ld2.remove(old_ld2[j])
                    #pass


            #pass
    return ld1,ld2
"""
def removeRedundants(chain_num,ld1): #### Conditions ####problem removing from ld1, ld2 but still existing in old
    #create a new list and instead of removing, add to list
    new_ld=[]
    cnL=[]
    i=0
    #print ("LD1_=",ld1)
    #print ("LEN_ld1=",len(ld1))
    while i<(len(ld1)-1):
        #print("LD1=",ld1[i],"i=",i)
        if ld1[i][0:4] == ld1[i+1][0:4]:
            cnL.append(ld1[i])
            #print ("NNNNNNNNNNNNNN",cnL)
            i+=1
            #pass
        else:
            cnL.append(ld1[i])
            #print ("NNNNNNNNNNNNNN",cnL)
            #print("LD22222=",ld1[i])
            if chain2num_dict[chain_num]<=len(cnL):
                new_ld.append(cnL[chain2num_dict[chain_num]-1])
            cnL=[]
            i+=1
    cnL.append(ld1[i])
            #print ("NNNNNNNNNNNNNN",cnL)
    #print("LD22222=",ld1[i])
    if chain2num_dict[chain_num]<=len(cnL):
        new_ld.append(cnL[chain2num_dict[chain_num]-1])

    return new_ld

def getLengthFromLengthDict(ld):
    l=[]
    for i in range(len(ld)):
        dct=ld[i]
        for key, val in dct.items():
            l.append(val)
    if max(l)!=min(l): 
        print ("Length Mismatch in Similar Alignments")
        sys.exit()
    firstName=list(ld[0].keys())[0]
    return firstName, max(l)

def findUniqueNodes(d_list): # browses through all the alignment indexes and creates a list of 
                             #sorted non-redundant nodes
    d_lst=[]
    dict_keys_dict={}
    #dict_keys=[]
    for d_ in d_list:
        #print(type(d_),"11111111111111",list(d_[0].keys())[0], d_[0])
        first_key=list(d_[0].keys())[0]
        dict_keys=[]
        #d_list.append(d_[])
        #print (d_[1])
        #print(type(d_[1]))
        for i in range(len(d_)):
            #print (d_[i], type(d_[i]),len(d_[i]))
            s=list(d_[i].keys())[0]
            #s=list(d_[i].keys())[0]
            #print(s, type(s))
            dict_keys.append(s)
            
            d_lst.append(s[0:4])
            #print(list(d_[i].keys())[0])
        dict_keys.remove(first_key)
        #dict_keys=removeRedundants(sorted(dict_keys))
        dict_keys_dict[first_key]=sorted(dict_keys)
    
    
    return list(set(d_lst)),dict_keys_dict

def convertMSA2PIR(MSAList):
    file_contents=[]
    print ("Preparing PIR contents ...")
    for i in range(len(MSAList)-1):
        row=MSAList[i]
        #print (row)
        #print(type(row))
        master_key=list(row.keys())[0]
        sub_list=list(row[master_key].keys())
        #(">P1;"+rname[0:4]+"\n")
        #print (master_key, sub_list, type(sub_list))
        print (">P1;"+master_key+"\n")
        file_contents.append(">P1;"+master_key+"\n")
        tag="structure:"+master_key+"_"
        sequence=""
        for chain_pdbs in sub_list:
            tag+=chain_pdbs[4]
            sequence+=row[master_key][chain_pdbs]
            if sequence.endswith("*"): sequence=sequence[:-1]+"/"
        if sequence.endswith("/"):sequence=sequence[:-1]+"*\n"
        tag+=": : : : : : : : \n"

        print (tag)
        print (sequence)
        file_contents.append(tag)
        file_contents.append(sequence)
        #print ("i=",i)

        #break
    row=MSAList[-1]
    master_key=list(row.keys())[0]
    sub_list=list(row[master_key].keys())
    print (">P1;"+master_key+"\n")
    file_contents.append(">P1;"+master_key+"\n")
    tag="sequence:"+master_key+"_"
    sequence=""
    print ("Chain_pdbs=",chain_pdbs)
    sys.exit()
    for chain_pdbs in sub_list:
            tag+=chain_pdbs[4]
            sequence+=row[master_key][chain_pdbs]
            if sequence.endswith("*"): sequence=sequence[:-1]+"/"
    if sequence.endswith("/"):sequence=sequence[:-1]+"*\n"
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
for file in inp_pir_list:
    name=os.path.basename(file).replace(".pir","")
    #print (tmpdir+name+"/"+name+"_aln.pir","77777777777777777")
    dd,LL=createDict(tmpdir+name+"/"+name+"_aln.pir")
    #print(dd,"88888888888888888")
    d.append(dd)
    fn,L_nm=getLengthFromLengthDict(LL)
    length_dict[fn]=L_nm
    #d.append(createDict(tmpdir+name+"/"+name+"_aln.pir"))
    #d.append(createDict(tmpdir+name+"/"+name+"_aln.pir"))
    #length_dict[name]=lengthOfAlignedFasta(d[-1])

print (length_dict)
print("########################################",name)

#sys.exit()

k_d_key=sorted(list(length_dict.keys()))
length_dict_list=[]
for k_d in k_d_key:
    length_dict_list.append(length_dict[k_d])

#sys.exit()

#dList is the list of the unique PDB codes in all the alignment files. 
#dict_keys_dict is the dictionary where the query sequence points to a list of the other alignment nodes

dList, dict_keys_dict=findUniqueNodes(d) 
print ("@@@@@@=",len(dList))
MAX_ALN_NUM=len(dList)
print (MAX_ALN_NUM)
print ("dList=",dList,"*****************************************\n")
print (dict_keys_dict)
print ("Chains=",fileNumber)
#print (dList[1])
#print (dList[2])
#print (dict_keys_dict)

print ("KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK")
#sorted keys is the query key of the chains
sorted_keys=sorted(list(dict_keys_dict.keys()))
print("111111")
print (sorted_keys)
print("111111")
for sorted_k in sorted_keys:
    dict_keys_dict[sorted_k]=sorted(dict_keys_dict[sorted_k])

#print ("sorted dict_key_dict=",dict_keys_dict)
#sys.exit()


#Create An all alignment dictionary for simplicity. Use the dict_keys_dict

aln_dict={}

for q_k_i in range(len(sorted_keys)):
    q_k=sorted_keys[q_k_i]
    aln_keys=dict_keys_dict[q_k]
    for _ii_ in range(len(d[q_k_i])):
        if q_k == list(d[q_k_i][_ii_].keys())[0]:
            aln_dict[q_k]=d[q_k_i][_ii_][q_k]
        for aln_k in aln_keys:
            if aln_k == list(d[q_k_i][_ii_].keys())[0]:
                aln_dict[aln_k]=d[q_k_i][_ii_][aln_k]
    


print("ALN_DICT=",aln_dict)
#sys.exit()


#The following is used to retrieve the chains.
query_k_list=[]
for query_k in dict_keys_dict:
    query_k_list.append(query_k)
    #print (query_k)
chain_dict[query_k[0:4]]=query_k_list
chain_dict_list.append({query_k[0:4]:chain_dict[query_k[0:4]]})
#print (sorted_keys)

dList_copy=list.copy(dList)
dList.remove(query_k[0:4])

print (chain_dict)

#sys.exit()


print(dList)
#sys.exit()




#The following is a dictionary of the 4 leter PDB codes mapping to the chains.
#Don't ask me how I did this.
"""
for i in range(len(dList)):
    index_k_list=[]
    #print ("HERE")
    print ("dList",i,dList[i])
    alreadyPassed=False
    for query_k in sorted_keys:
        chained_pdbs_list=dict_keys_dict[query_k]
        #print ("_CHAINED_PDBS=",chained_pdbs)
        #index_k_list.append(chained_pdbs)
        for chained_pdbs in chained_pdbs_list:
            #if chained_pdbs[0:4] == '6MWC': print ("FLAW=",chained_pdbs,":")
            if dList[i]==chained_pdbs[0:4] and alreadyPassed==False:
                print ("HERE")
                for _i in range (fileNumber):
                    _chain=chained_pdbs[4]
                    if chained_pdbs[0:4] == '6MWC': print ("_i value =",_i,"FLAW=",chained_pdbs,":",_chain)
                    chainedname=chained_pdbs[0:4]+reverse_chain2num_dict[chain2num_dict[_chain]+_i]
                    chainedname=chained_pdbs
                    print ("_i=",_i,"_CHAIN=",_chain,"Sum=",chain2num_dict[_chain]+_i,"Chainedname=",chainedname)
                    index_k_list.append(chainedname)
                chain_dict[dList[i]]=index_k_list
                chain_dict_list.append({dList[i]:chain_dict[dList[i]]})
                alreadyPassed=True
                break
            else:
                continue
    #print ("i==",i,sorted_keys[i])
    #chain_dict[dList[i]]=index_k_list
"""

#print (sorted_keys)
#sys.exit()

# The following is an attempt on a better version:

new_dict_keys_dict={}

for i in range(len(dList)):
    index_k_list=[]
    #print ("HERE")
    print ("dList",i,dList[i])
    alreadyPassed=False
    for query_k in sorted_keys:
        chained_pdbs_list=sorted(dict_keys_dict[query_k])
        #print ("_CHAINED_PDBS=",chained_pdbs)
        #index_k_list.append(chained_pdbs)
        chain_count=0
        print ("chained_pdbs_list=",chained_pdbs_list, len(chained_pdbs_list))
        new_chained_pdb_list=[]
        for chained_pdbs in chained_pdbs_list:
            print ("chained_pdbs=",chained_pdbs)
            #sys.exit()
            #if chained_pdbs[0:4] == '6MWC': print ("FLAW=",chained_pdbs,":")
            if dList[i]==chained_pdbs[0:4] and chain_count<fileNumber: # and alreadyPassed==False:
                    print ("HERE")
                    new_chained_pdb_list.append(chained_pdbs)
                    _chain=chained_pdbs[4]
                    print (chained_pdbs)
                    #new_dict_keys_dict[]
                    if chained_pdbs[0:4] == '6MWC': print ("_i value =","FLAW=",chained_pdbs,":",_chain)
                    #chainedname=chained_pdbs[0:4]+reverse_chain2num_dict[chain2num_dict[_chain]+_i]
                    chainedname=chained_pdbs
                    #print ("_i=",_i,"_CHAIN=",_chain,"Sum=",chain2num_dict[_chain]+_i,"Chainedname=",chainedname)
                    index_k_list.append(chainedname)
                    chain_dict[dList[i]]=new_chained_pdb_list #index_k_list
                    #chain_dict_list.append({dList[i]:chain_dict[dList[i]]})
                    chain_dict_list.append({dList[i]:chain_dict[dList[i]]})
                    alreadyPassed=True
                    chain_count+=1
                #break
            else:
                #break
                continue
print ("chained_dict",chain_dict)
initializing_chain_dict={}





sys.exit()    

print (chain_dict)
print (len(chain_dict))
print (len(chain_dict_list))
firstOnList=chain_dict_list[0]
chain_dict_list.remove(firstOnList)
chain_dict_list.insert(len(chain_dict_list),firstOnList)
print (chain_dict_list)
print (len(chain_dict_list))
#sys.exit()


# Generate MSA graph and initialize nodes with gaps. Use the dList
print ("Generate MSA graph and initialize nodes with gaps")
print (length_dict_list)
MSAGraphMatrix_dict=[]
MSAGraphMatrix_list=[]
for row in chain_dict_list:
    print ("Row=",row)
    for key, value in row.items():
        print ("K,V=",key,value)
        l_msa=[]
        l_msa_dict={}
        for chains in range(len(value)):
            l_msa.append({value[chains]:"-"*(length_dict_list[chains]-1)})
            l_msa_dict[value[chains]]="-"*(length_dict_list[chains]-1)
        MSAGraphMatrix_dict.append({key:l_msa})
        MSAGraphMatrix_list.append({key:l_msa_dict})
    #break
#sys.exit()
#print (MSAGraphMatrix_dict)
#print (MSAGraphMatrix_list)
#sys.exit()

#print (d)
#print (len(d[2]))

#print (type(d[0][1]))
#print(len(d[0]))
#print (d[0])
#print (d[0][0])
#sys.exit()



# The following parses through each alignment file and removes redundant chains 
# or homomeric chains that belong to other subunits

for i in range(len(sorted_keys)):
    d_1=dict_keys_dict[sorted_keys[i]]
    #d_2=dict_keys_dict[sorted_keys[i+1]]
    print(sorted_keys[i][-1])
    #sys.exit()
    new_d_1=removeRedundants(sorted_keys[i][-1],d_1)
    dict_keys_dict[sorted_keys[i]]=new_d_1
    #print (dict_keys_dict)
    #sys.exit()
print (dict_keys_dict)

#Create An all alignment dictionary for simplicity. Use the dict_keys_dict

aln_dict={}

for q_k_i in range(len(sorted_keys)):
    q_k=sorted_keys[q_k_i]
    aln_keys=dict_keys_dict[q_k]
    for _ii_ in range(len(d[q_k_i])):
        if q_k == list(d[q_k_i][_ii_].keys())[0]:
            aln_dict[q_k]=d[q_k_i][_ii_][q_k]
        for aln_k in aln_keys:
            if aln_k == list(d[q_k_i][_ii_].keys())[0]:
                aln_dict[aln_k]=d[q_k_i][_ii_][aln_k]
    


#print(aln_dict)


#Update the MSAs to generate the final concatenated MSA

#print(d)
#print (type(d))
#print(len(d))
#print (MSAGraphMatrix_list)
#print (MSAGraphMatrix_list[0])
#print(len(MSAGraphMatrix_list))
#print (chain_dict)
#print (type(chain_dict))
#sys.exit()

print ("$$$$START####")
# Read the chain_dict. And update the MSAGraphMatrix_list accordingly
toUpdate_key_list=[]
master_keys_list=list(chain_dict.keys())
for master_key in master_keys_list:
    #print (master_key)
    child_key_list=chain_dict[master_key]
    print (child_key_list)
    toUpdate_key_list+=child_key_list
    pass
print(toUpdate_key_list)
#sys.exit()
print (MSAGraphMatrix_dict)
print (len(MSAGraphMatrix_dict))
print (type(MSAGraphMatrix_dict))

print (MSAGraphMatrix_list)
print (len(MSAGraphMatrix_list))
print (type(MSAGraphMatrix_list))
#for update_key in toUpdate_key_list:
updatedMSAList=[]
print ("****************")
for update_key in toUpdate_key_list:
    for _ii in range(len(MSAGraphMatrix_list)):
        row = MSAGraphMatrix_list[_ii]
        #print (row[list(row.keys())[0]])
        #print (type(row[list(row.keys())[0]]))
        #print (len(row['5GNI']))
        #break
        #row = MSAGraphMatrix_list[_ii]
        if update_key[0:4] == list(row.keys())[0]:
            if aln_dict.get(update_key)!=None:
                MSAGraphMatrix_list[_ii][update_key[0:4]][update_key]=aln_dict[update_key]
                break

print(MSAGraphMatrix_list)
file_contents=convertMSA2PIR(MSAGraphMatrix_list)


#print (outdir+"pir_out/"+name)
#sys.exit()
if not os.path.isdir(outdir+"pir_out/"): os.makedirs(outdir+"pir_out/")
with open (outdir+"pir_out/"+name[0:4]+"_0.pir","w") as f:
    f.writelines(file_contents)
sys.exit()


"""
for _kk in range(len(MSAGraphMatrix_list)):
    row=MSAGraphMatrix_list[_kk]
    print (_kk,"Row==",row)
    #print("Type_row=",type(row)) #dict
    #print("L_row=",len(row))    #chains
    main_pdb_key=list(row.keys())[0]
    print ("Main_key=",main_pdb_key)
    print("Type_main_key=",type(row[main_pdb_key]))
    print("L_main_key=",len(row[main_pdb_key]))

    sub_key_list=sorted(list(row[main_pdb_key].keys()))
    print ("sub_key_list=",sub_key_list)
    #break
    for _ii in range(0,len(sub_key_list)):
        each_key=sub_key_list[_ii]
        #print (each_key,d[_ii][0].keys())
        print (main_pdb_key, each_key)
        print (d[_ii][0])
        #break
        #print (d[_ii][0][each_key])
        if d[_ii][0].get(each_key) != None:
            MSAGraphMatrix_list[_kk][main_pdb_key][each_key]=d[_ii][0][each_key]
            #print (d[_ii][0][each_key])


    break
"""
"""
for _kk in range(len(MSAGraphMatrix_list)):
    row=MSAGraphMatrix_list[_kk]
    print (_kk,"Row==",row)
    #print("Type_row=",type(row)) #dict
    #print("L_row=",len(row))    #chains
    main_pdb_key=list(row.keys())[0]
    print ("Main_key=",main_pdb_key)
    print("Type_main_key=",type(row[main_pdb_key]))
    print("L_main_key=",len(row[main_pdb_key]))

    sub_key_list=sorted(list(row[main_pdb_key].keys()))
    print ("sub_key_list=",sub_key_list)
    #break
    for _ii in range(0,len(sub_key_list)):
        each_key=sub_key_list[_ii]
        #print (each_key,d[_ii][0].keys())
        print (main_pdb_key, each_key)
        print (d[_ii][0])
        #break
        #print (d[_ii][0][each_key])
        if d[_ii][0].get(each_key) != None:
            MSAGraphMatrix_list[_kk][main_pdb_key][each_key]=d[_ii][0][each_key]
            #print (d[_ii][0][each_key])


    break
"""
#print (MSAGraphMatrix_list)
#for chain_stack_i in d:
#    print (chain_stack_i)


sys.exit()
#print ("UUUUUUUUUUUUUUUUUUUUUUUUUUUUUPPPPPPPPDDDDAAAAAATTTTTTTTEEEEEE")



#print (dict_keys_dict)
"""
chain_dict=generateChainDict(d)

combo_d={}

keys=list(d[0].keys()) #sorted
for i in range(1,len(d)-1):
    for k in keys:
        d_k=findInterface(k,d[i])
        joinInterface(k,d_k)

"""



# Separate script. Pass individual pir files to the 
# Use Modeller to generate rolling templates



