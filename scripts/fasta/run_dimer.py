import os, sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

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

fasta_file_A=os.path.abspath("/data/farhan/multicom4s/tool/multicom4s_TMB/scripts/homotrimer/homotrimer_fasta/1RERA.fasta")
fasta_file_B=os.path.abspath("/data/farhan/multicom4s/tool/multicom4s_TMB/scripts/homotrimer/homotrimer_fasta/1RERB.fasta")
output_folder=os.path.abspath("./out")+"/"
template_fasta_folder=os.path.abspath("./")
if not os.path.isdir(output_folder): os.makedirs(output_folder)

local_aln=[]
orig_aln=[]

fasta_A=readFastaFile(fasta_file_A)
fasta_B=readFastaFile(fasta_file_B)
print (os.path.basename(fasta_file_A))
name=os.path.basename(fasta_file_A).replace(".fasta","")+"_"+os.path.basename(fasta_file_B).replace(".fasta","")

rname="1LD4A_1LD4B"

template_name_A=rname.split("_")[0]
template_name_B=rname.split("_")[1]
fasta_template_A=readFastaFile(template_fasta_folder+"/"+template_name_A+".fasta")
fasta_template_B=readFastaFile(template_fasta_folder+"/"+template_name_B+".fasta")
print ("Orig Fastas A")
print (fasta_A)
print ("Template Fastas A")
print (fasta_template_A)
print ("Orig Fastas B")
print (fasta_B)
print ("Template Fastas B")
print (fasta_template_B)
alignments_A=pairwise2.align.globalms(fasta_A, fasta_template_A,5,-4,-1,-0.1)
alignments_B=pairwise2.align.globalms(fasta_B, fasta_template_B,5,-4,-1,-0.1)
best_aln_A=selectBestAlignment(alignments_A)
best_aln_B=selectBestAlignment(alignments_B)

print (best_aln_A[0])
print (best_aln_A[1])
print ("########################3")

print (best_aln_B[0])
print (best_aln_B[1])
orig_aln.append(best_aln_A[0]+"/"+best_aln_B[0]+"*")
local_aln.append(best_aln_A[1]+"/"+best_aln_B[1]+"*")

print (len(local_aln),len(orig_aln))
#sys.exit()
with open (output_folder+"/"+name+"_"+str(0)+".pir","w") as f:
        f.write(">P1;"+rname[0:4]+"\n")
        chain_A=rname.strip().split("_")[0][-1]
        chain_B=rname.strip().split("_")[1][-1]
        f.write("structureN:"+rname[0:4]+"::"+chain_A+"::"+chain_B+"::::\n")
#        f.write(local_aln[i]+"\n")
        f.write(local_aln[0]+"\n")
        f.write("\n")
        f.write(">P1;"+name+"\n")
        targ_chain_A=name.strip().split("_")[0][-1]
        targ_chain_B=name.strip().split("_")[1][-1]
        f.write("sequence: :"+targ_chain_A+" : :"+targ_chain_B+" : : : : : \n")
#        f.write(orig_aln[i]+"\n")
        f.write(orig_aln[0]+"\n")





