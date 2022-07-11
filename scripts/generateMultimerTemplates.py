#this script takes input the fasta files of the templates according to the chains
#usage: python generateMultimerTemplates.py <input_fasta_directory> <output_dir> <--hhblits>
#1. Generates monomer templates-
#      (a) By aligning using hhblits, bfd and a custom monomer hhblits searchable database.
#      (b) Using the PIR alignment database.
#
#2. Generates an alignment file of all the templates in the PIR format. The alignments are separated by '/'
#3. Runs modeller on the file in (2)


import os, sys
from glob import glob

class Node:
    def __init__(self,pir_indx,pdb_code,chained_pdb_code,L,points2, prev1, target):
        self.pir_idx=pir_idx
        self.pdb=pdb_code
        self.chain_pdb=chained_pdb_code
        self.chain=chain_pdb[4]
        self.length=L
        self.partner2=points2
        self.partner_prev=prev1
        self.target=target
        self.aln="-" * L
    
    def __str__(self):
        s="Original PDB Code: "+self.code()+"\n"
        ch="Chain is: "+self.chain()+"\n"
        ln="Length is: "+self.length()+"\n"
        p2="Next Partner: "+self.partner2()+"\n"
        p_prev="Previous Partner: "+self.partner_prev()+"\n"
        isT="Is Target: "+self.target()+"\n"
        al="Alignment:\n"+self.alignment()+"\n"
        return s+ch+ln+p2+p_prev+isT+al

    def pir_idx(self):
        return self.pir_idx
    def code(self):
        return self.pdb
    def chain(self):
        return self.chain
    def chained_pdb(self):
        return self.chain_pdb
    def target(self):
        return self.target
    def length(self):
        return self.length
    def partner2(self):
        return self.partner2
    def partner_prev(self):
        return self.partner_prev
    def alignment(self):
        return self.aln
    def replaceAln(self,alignment):
        self.aln=alignment






length_dict={}
chain_dict={}
s='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
chain2num_dict={}
reverse_chain2num_dict={}
#global MAX_ALN_NUM
MAX_ALN_NUM=0
for i in range(len(s)):
    chain2num_dict[s[i]]=i+1
    reverse_chain2num_dict[i+1]=s[i]




def prob(e):
    return e[2]

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
    print (d)
    
    m=0
    for k in d:
        #keys=k.keys()
        m=len(k)
        if MAX_ALN_NUM<m: MAX_ALN_NUM=m
    """
    for k in d:
        keys=k.keys()
        m=len(keys)
        if MAX_LENGTH<m: MAX_LENGTH=m
    """
    return d, L_chains
                

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

def generateChainDict(d):
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
def removeRedundants(ld1,ld2): #### Conditions ####problem removing from ld1, ld2 but still existing in old
    #create a new list and instead of removing, add to list
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

def findUniqueNodes(d_list):
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
            print (d_[i], type(d_[i]),len(d_[i]))
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


#########################################################################################3
# Main
##########################################################################################
inpdir=os.path.abspath(sys.argv[1])+"/"
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



fastalist=sorted(glob(inpdir+"*.fasta"))
print (fastalist)
for file in fastalist:
    pir_name=fastapirdir+os.path.basename(file).replace(".fasta",".pir")
    print (pir_name)
    inp_pir_list.append(pir_name)
    os.system("python fasta2pir.py "+file+" "+pir_name)

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

d=[]
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
print("########################################")
dList, dict_keys_dict=findUniqueNodes(d) 
print (len(dList))
print (dList,"*****************************************")
#print (dList[1])
#print (dList[2])
print (dict_keys_dict)
sorted_keys=sorted(list(dict_keys_dict.keys()))
print (sorted_keys)

for i in range(len(sorted_keys)-1):
    d_1=dict_keys_dict[sorted_keys[i]]
    d_2=dict_keys_dict[sorted_keys[i+1]]
    new_d_1,new_d_2=removeRedundants(d_1,d_2)
    dict_keys_dict[sorted_keys[i]]=new_d_1
    dict_keys_dict[sorted_keys[i+1]]=new_d_2

print ("UUUUUUUUUUUUUUUUUUUUUUUUUUUUUPPPPPPPPDDDDAAAAAATTTTTTTTEEEEEE")

print (dict_keys_dict)
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



