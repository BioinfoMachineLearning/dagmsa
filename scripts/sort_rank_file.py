#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:01:04 2019

@author: farhan
"""
#sorts an template rank file and prints its contents
#usage: python sort_rank_file.py <rank_file.hhr> <reverse=True> 
import os,sys

def prob(e):
    return e[-1]


rankfile=os.path.abspath(sys.argv[1])
rev=sys.argv[2]
label=""

def readFile(file):
    contents=[]
    fasta=""
    line=""
    lab=""
    with open (file,"r") as f:
        for line in f:
            if (line.strip().startswith("Ranked")):
                lab=line
                continue
            split=line.strip().split()
            #contents.append([split[0],split[1],split[2],float(split[3])])
            contents.append([str(split[1]+"\t"+split[2]),float(split[3])])
            
    #if (line.strip()=="END"): fasta=fasta[0:len(fasta)-3]
    return lab,contents

label,contents=readFile(rankfile)
def toNum(contents):
    new_contents=[]
    for line in contents:
        split=line.split()
        i=int(split[0])
        j=int(split[1])
        zero=int(split[2])
        dist=int(split[3])
        val=float(split[4])
        #print (i,j,zero,dist,val)
        new_contents.append([i,j,zero,dist,val])
    return new_contents



#contents=toNum(contents)
if rev.strip()=="True": contents.sort(key=prob,reverse=True)
if rev.strip()=="False": contents.sort(key=prob,reverse=False)

#print (contents)
#contents.insert(0,label)
print(label)
#print(fasta)
i=0
for item in contents:
    i+=1
    print(str(i)+"\t"+str(item[0])+"\t"+str(item[1])+"\n")#+" "+str(item[5]))
#print("END")

