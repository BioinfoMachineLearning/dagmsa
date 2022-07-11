# Comparative modeling by the AutoModel class
#
# Demonstrates how to build multi-chain models, and symmetry restraints
#
#from ctypes.wintypes import tagRECT
#import tarfile
#from modeller import *
#from modeller.automodel import *    # Load the AutoModel class
import os, sys

#from unicodedata import name

#Pipeline:
#1. Read the main PIR file
#2. Create a new pir file containing the first two chains only. 
#3. Pass this to modeller in the script run_modeller.py
#4. Use #2 to combine the alignments 
#5. Combine the output PDB file from #3 to single chain and reindex.
#6. Rerun using same [recursive]


def merge_pair_of_aln_2_pir(first, second):
    chain_A=""
    chain_B=""
    file_contents=[]
    for i in range(len(first)-1):
        key_A=list(first[i][0].keys())[0]
        key_B=list(second[i][0].keys())[0]
        name=key_A[0:4]
        chain_A=key_A[4]
        chain_B=key_B[4]
        #print (">P1;"+name+"\n")
        file_contents.append(">P1;"+name+"\n")
        tag="structure:"+name+":"+chain_A+" : :"+chain_B+" : : : : : \n"
        sequence=first[i][0][key_A]+"/"+second[i][0][key_B]
        sequence=sequence.replace("*","")+"*\n\n"
        file_contents.append(tag)
        file_contents.append(sequence)
    
    key_A=list(first[-1][0].keys())[0]
    key_B=list(second[-1][0].keys())[0]
    name=key_A[0:4]
    chain_A=key_A[4]
    chain_B=key_B[4]
    #print (">P1;"+name+"\n")
    file_contents.append(">P1;"+name+"\n")
    tag="seequence:"+name+":"+chain_A+" : :"+chain_B+" : : : : : \n"
    sequence=first[-1][0][key_A]+"/"+second[-1][0][key_B]
    sequence=sequence.replace("*","")+"*\n\n"
    file_contents.append(tag)
    file_contents.append(sequence)

    return file_contents

def prep_paired_pir(aln_dict_list,aln_list,pdb_files_folder):
    max_chain=len(aln_dict_list)
    temp_pdb_files_folder=outdir+"template_atom_files/"
    if not os.path.isdir(temp_pdb_files_folder): os.makedirs(temp_pdb_files_folder)
    files_list=list(aln_list[0].keys())
    print (files_list)
    #sys.exit()
    not_found_count=0
    for file in files_list:
        #exit_code=os.system("scp "+pdb_files_folder+file[0:4]+"* "+temp_pdb_files_folder)
        print ("HERE")
        exit_code=os.system("python createPDBFolder_v2.py "+file[0:4]+" "+pdb_files_folder+" "+temp_pdb_files_folder)
        #sys.exit()
        if exit_code!=0:
            not_found_count+=1
            #exit_code_2=os.system("python createPDBFolder_v2.py "+file[0:4]+" "+pdb_files_folder+" "+temp_pdb_files_folder)
            #if exit_code_2!=0:
            #    not_found_count+=1
        

    print ("Not_found_count=",not_found_count)

    print ("MAX_CHAIN=",max_chain)
    for i in range(max_chain-1):
        first_chain=aln_dict_list[i]
        sec_chain=aln_dict_list[i+1]
        #print (first_chain[-1])
        #print ("Type=",type(first_chain),"Len=",len(first_chain))
        #files_list_B=list(aln_list[i+1].keys())
        file_contents=merge_pair_of_aln_2_pir(first_chain,sec_chain)
        #print (type(first_chain[-1][0]))
        first_name=list(first_chain[-1][0].keys())[0][0:4]
        first_chain_name=list(first_chain[-1][0].keys())[0][4]
        second_chain_name=list(sec_chain[-1][0].keys())[0][4]
        #print("CHAIN_NAMES=",first_chain_name,second_chain_name)
        #sys.exit()
        temp_aln_pir_file=outdir+first_name+"_"+first_chain_name+second_chain_name+".temp.pir"

        with open (temp_aln_pir_file,"w") as f:
            f.writelines(file_contents)
        
        #Run Modeller
        #os.system("python run_modeller_v2.py "+temp_aln_pir_file+" "+temp_pdb_files_folder)
        
        break



        
    return


chain_number=0


# directories for input atom files
aln_pir_file=os.path.abspath(sys.argv[1])

pdb_files_folder=os.path.abspath(sys.argv[2])+"/"

chain_number=int(sys.argv[3])

outdir=os.path.abspath(sys.argv[4])+"/"

if not os.path.isdir(outdir): os.makedirs(outdir)

alignment_list=[]
alignment_dict_list=[]

name=""

first_aln_dict={}
first_aln_index_list=[]



for i in range(chain_number):
    alignment_list.append({})
    alignment_dict_list.append([])

print (alignment_list)



knowns_list=[]
kns="" #tuple for multiple templates
query=""
with open (aln_pir_file) as f:
    for line in f:
        if line.startswith(">"):
            line=f.readline().strip()
            split=line.strip().split(":")[1]
            #kns=split[0:4]
            pdb_names=split.strip().split("_")[0] #.upper()
            #pdb_names=split.strip().split("_")[0] #normal
            chain_series=split.strip().split("_")[1]
            aln=f.readline().strip().split("/")
            #print (len(alignment_list),len(chain_series),chain_series)
            for _i in range(len(chain_series)):
                #print ("_i=",_i)
                knowns_list.append(pdb_names+chain_series[_i])
                alignment_list[_i][pdb_names+chain_series[_i]]=aln[_i]
                alignment_dict_list[_i].append([{pdb_names+chain_series[_i]:aln[_i]}])

                #break
            #continue

#print (alignment_list)
print (len(alignment_list))

#print(len(alignment_list[-1]))
#print (alignment_list[0])

print (len(alignment_dict_list[0]))
#print (alignment_dict_list[0])

prep_paired_pir(alignment_dict_list,alignment_list,pdb_files_folder)

sys.exit()

local=pdb_files_folder

if len(knowns_list) == 2:
    query=knowns_list[-1]
    kns=knowns_list[0]
elif len(knowns_list)>2:
    query=knowns_list[-1]
    copylist=knowns_list.copy()
    copylist.remove(query)
    kns=tuple(copylist)
else:
    sys.exit("No templates found in list")

print (query)
print (kns)
#sys.exit()

#env.io.atom_files_directory = ['.', pdb_files_folder]
env.io.atom_files_directory = ['',local]
print (env.io.atom_files_directory)

# Be sure to use 'MyModel' rather than 'AutoModel' here!
a = MyModel(env,
            alnfile  = aln_pir_file ,     # alignment filename
            knowns   = kns,              # codes of the templates
            sequence = query
            )              # code of the target

a.starting_model= 1                # index of the first model
a.ending_model  = 1                # index of the last model
                                   # (determines how many models to calculate)
a.make()                           # do comparative modeling