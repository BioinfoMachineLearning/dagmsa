import os, sys
from glob import glob
from OneMSASubGraph import OneMSASubGraph

class DAGMSA:

    def __init__(self, onemsa_list):
        #One MSA list files. Convert to OneMSASubGraph format
        self.onemsa_list=[]
        for onemsa in onemsa_list:
            self.onemsa_list.append([onemsa, onemsa.getAvgEvalue])
        pass
