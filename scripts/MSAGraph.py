
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


class SubGraph:

    def __init__(self, number, pdb_code, node_list):
        pass

class MSAGraph:
    def __init__(self, sub_graph_list):
        pass