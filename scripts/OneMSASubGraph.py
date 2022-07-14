import os, sys
from glob import glob
from Node import Node
class OneMSASubGraph:

    def __init__(self, number, pdb_code, node_list):
        # we need to generate the nodes and intialize from the list. How to do it? ???????
        stringofchains="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        stringofchains+=stringofchains.lower()
        stringofchains+="123456789"
        self.node_list=node_list
        self.avgevalue=0
        for node in node_list:
            self.avgevalue+=node.getEvalue()
        if len(node_list)!=0: self.avgevalue/=len(node_list)

        pass

    def getAvgEvalue(self):
        return self.avgevalue

    def getAlignment(self,format="PIR",separator="/"):
        if format.upper()=="PIR":
            pass
        elif format.upper()=="AF":
            pass

        pass