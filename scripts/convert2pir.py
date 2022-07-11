#this script takes the rank file, local file and the fasta files as input and generates a .pir file for modeller
import os, sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

rank_file=os.path.abspath(sys.argv[1])
local_file=os.path.abspath(sys.argv[2])
fasta_file_A=os.path.abspath(sys.argv[3])
fasta_file_B=os.path.abspath(sys.argv[4])
output_folder=os.path.abspath(sys.argv[5])

def readRankFile(file):
    rank=[]
    i=0
    with open (file) as f:
        for line in f:
            if i==10: break;
            if line.startswith("Ranked"): continue
            split=line.strip().split()
            rank.append([split[1],split[2],split[3]])
            i+=1
    return rank

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

rank=readRankFile(rank_file)
#print (len(rank))
#print (rank[0])
#print (rank[9])

fasta_A=readFastaFile(fasta_file_A)
fasta_B=readFastaFile(fasta_file_B)
print (os.path.basename(fasta_file_A))
name=os.path.basename(fasta_file_A).replace(".fasta","")+"_"+os.path.basename(fasta_file_B).replace(".fasta","")
rank_list=[]
for rr in rank:
    rank_list.append(rr[0])

print (len(fasta_A)+len(fasta_B))
#print (len(fasta_B))

local_aln, orig_aln=readLocalFile(local_file,name,fasta_A+fasta_B, rank_list)


#print (fasta_A+"/"+fasta_B+"*\n")
for rname in rank_list:
    print (local_aln[rname])
    #local_aln[rname]=local_aln[rname][0:len(fasta_A)]+"/"+local_aln[rname][len(fasta_A):]+"*\n"
    #local_aln[rname]=local_aln[rname][0:len(fasta_A)]+local_aln[rname][len(fasta_A):]+"*\n"
    local_aln[rname]=local_aln[rname]+"\n"
    orig_aln[rname]=orig_aln[rname]+"\n"
    
    #print ("ORI:",orig_aln[rname])
    #print ("LOC:",local_aln[rname])
    


print (len(fasta_A+"/"+fasta_B+"*\n"))
for i in range(10):
    print (len(local_aln[rank_list[i]]))
    print (len(orig_aln[rank_list[i]]))
    print ("//" in local_aln[rank_list[i]])

for i in range(len(rank_list)):
    with open (output_folder+"/"+name+"_"+str(i)+".pir","w") as f:
        rname=rank_list[i]
        f.write(">P1;"+rname+"\n")
        f.write("sequence:"+rname+"::::::::\n")
        f.write(local_aln[rname])
        f.write("\n")
        f.write(">P1;"+name+"\n")
        f.write(" : : : : : : : : : \n")
        f.write(orig_aln[rname])








