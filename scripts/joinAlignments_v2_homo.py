#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 03:13:59 2020

@author: farhan
"""

#this script reads the two heteromeric alignments in .a3m format
#creates respectie dictionaries and indices
#Matches the Headers in the a3m files for interacting pairs
#then creates the concatenated alignment file in .a3m format

import os, sys, shutil
import numpy as np

def getName(file):
    return os.path.basename(file.strip()).split(".")[0]

def createDictFromA3m(file,name):
    aln_dict={}
    with open (file) as f:
        for line in f:
            if line.startswith(">"):
                idx=line.strip().replace(">","")
                #idx=idx.split("/")[0]
                if name in idx: idx=name
                line=f.readline()
                aln_dict[idx]=line.strip()            
    return aln_dict

def saveFastaDictionary(filename,dict):
    with open (filename,"w") as f:
        for k,v in dict.items():
            f.write(k+":"+v+"\n")

def loadFastaDictionary(filename):
    d={}
    with open (filename) as f:
        for line in f:
            idx=line.strip().split(":")[0].strip()
            fasta=line.strip().split(":")[1].strip()
            d[idx]=fasta
    return d

def getPPIDictFromPPIFile(_file):
   d = {}
   with open (_file) as f:
      for line in f:
          key=line.strip().split(":")[0]
          sp=line.strip().split(":")[1].split(",")
          for i in range(len(sp)):
             sp[i]=sp[i].strip()
          d[key]=sp
   return d



#aln_A=
#aln_B=
aln_A=os.path.abspath(sys.argv[1])
aln_B=os.path.abspath(sys.argv[2])
#ppi_file=os.path.abspath(sys.argv[3])
outdir=os.path.abspath(sys.argv[3])+"/"
dbtypeflag=sys.argv[4]

if dbtypeflag!="-string" or dbtypeflag!="-bfd" or dbtypeflag!="-uniprot" or dbtypeflag!="-uniref":
    sys.exit("Please specify what type of IDs are used:-string,-bfd,-uniprot,-uniref")


missing_dir=outdir+"missing/"
print ("Output will be stored in "+outdir)
if not os.path.isdir(outdir): os.makedirs(outdir)
if not os.path.isdir(missing_dir): os.makedirs(missing_dir)
print ("Aln_file_A: "+aln_A)
print ("Aln_file_B: "+aln_B)
#loading the ppi_dict can fail. This file is very big. So loading it only once is logical

ppi_file = "/exports/store2/deepcomplex/datasets/ppi_dict_split00.txt" #np.load(ppi_file, allow_pickle='TRUE')

#ppi_dict = getPPIDictFromPPIFile(ppi_file)

#print (type(ppi_dict))
#print (ppi_dict[0])
#print (ppi_dict[1])
#print (ppi_dict[2])
#print (ppi_dict[-1])


#print ("ppi dict loaded!")
#sys.exit()
#failure=[]

if not os.path.exists(aln_A):
   with open ("failed_join.txt","a+") as fail:
       fail.write(aln_A+"\n")
   sys.exit(aln_A+" not found! Quitting")
if not os.path.exists(aln_B):
   with open ("failed_join.txt","a+") as fail:
       fail.write(aln_B+"\n")
   sys.exit(aln_B+" not found! Quitting")

#extract the names
name_A=getName(aln_A)
name_B=getName(aln_B)
#create the fasta dictionaries from the .a3m files
dict_A=createDictFromA3m(aln_A,name_A)
dict_B=createDictFromA3m(aln_B,name_B)
#temporary working directory
tmpdir=outdir+"tmpdir"+name_A+"_"+name_B+"/"
print (tmpdir)
#save the dictionaries as .txt files
if not os.path.isdir(tmpdir): os.makedirs(tmpdir)
dictfile_A=tmpdir+name_A+"_dict.txt"
dictfile_B=tmpdir+name_B+"_dict.txt"
saveFastaDictionary(dictfile_A,dict_A)
saveFastaDictionary(dictfile_B,dict_B)

keys_A=list(dict_A.keys())
keys_B=list(dict_B.keys())
keys_A.remove(name_A)
keys_B.remove(name_B)
pairs_AB=[]
pairs_AB.append([name_A,name_B])
#here comes the big search. This is the longest step. Needs to be optimized.
print ("Len (A)",len(keys_A))
print ("Len (B)",len(keys_B))

        #here comes the big search. This is the longest step. Needs to be optimized.
if (True):
        if os.path.exists(missing_dir+name_A+"_"+name_B+"_missing.txt"): os.remove(missing_dir+name_A+"_"+name_B+"_missing.txt")
        for idx_A in keys_A:
            interactions_with_idx_A=[]
            print ("Processing idx_A: "+idx_A)
            #get the interacting proteins for this idx
            if (idx_A.startswith(">tr|")):
                pass
            elif (idx_A.startswith(">Uniref") or idx_A.startswith(">UniRef")):
                pass
            elif (idx_A.startswith(">Uniref")):
                pass
            for idx_B in keys_B:
                if idx_A.split("/")[0] == idx_B.split("/")[0]:
                    pairs_AB.append([idx_A,idx_B])
            

        if os.path.exists(tmpdir+"pairs_AB.a3m"): os.remove(tmpdir+"pairs_AB.a3m")

        print ("Total of ",len(pairs_AB),"alignments found!")

        with open (tmpdir+"pairs_AB.txt","w") as f:
            for pair in pairs_AB:
                f.write(pair[0]+","+pair[1]+"\n")

        with open (tmpdir+"pairs_AB.a3m","w") as fa3m:
            for pair in pairs_AB:
                fa3m.write(">"+pair[0]+"_"+pair[1]+"\n")
                fa3m.write(dict_A[pair[0]].strip()+dict_B[pair[1]].strip()+"\n")

#        shutil.copy2(tmpdir+"pairs_AB.a3m",outdir+name_A+"_"+name_B+".a3m")

#        os.system("grep -v '^>' "+outdir+name_A+"_"+name_B+".a3m"+" | sed 's/[a-z]//g' > "+outdir+name_A+"_"+name_B+".aln")

#os.system("rm -rf "+tmpdir)














