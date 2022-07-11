# Comparative modeling by the AutoModel class
#
# Demonstrates how to build multi-chain models, and symmetry restraints
#
#from ctypes.wintypes import tagRECT
#import tarfile
from modeller import *
from modeller.automodel import *    # Load the AutoModel class
import os, sys

log.verbose()

# Override the 'special_restraints' and 'user_after_single_model' methods:
class MyModel(AutoModel):
#    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
#       self.rename_segments(segment_ids=[’A’, ’B’],renumber_residues=[1, 1])

    def special_restraints(self, aln):
        # Constrain the A and B chains to be identical (but only restrain
        # the C-alpha atoms, to reduce the number of interatomic distances
        # that need to be calculated):
        s1 = Selection(self.chains['A']).only_atom_types('CA')
        s2 = Selection(self.chains['B']).only_atom_types('CA')
        #self.restraints.symmetry.append(Symmetry(s1, s2, 1.0))
    def user_after_single_model(self):
        # Report on symmetry violations greater than 1A after building
        # each model:
        self.restraints.symmetry.report(1.0)

env = Environ()
# directories for input atom files
aln_pir_file=os.path.abspath(sys.argv[1])

pdb_files_folder=os.path.abspath(sys.argv[2])

knowns_list=[]
kns="" #tuple for multiple templates
query=""
with open (aln_pir_file) as f:
    for line in f:
        if line.startswith(">"):
            split=line.strip().split(";")[1]
            #kns=split[0:4]
            knowns_list.append(split)
            #continue
#local="/data/farhan/SoftwareTools/DeepComplex/DeepComplex/scripts/database_generating_scripts/hetero_std/dimer_pdbs"
#local="/data/farhan/SoftwareTools/DeepComplex/DeepComplex/scripts/database_generating_scripts/hetero_std/hetero_atom_files_chain_added"
#local="/data/farhan/SoftwareTools/DeepComplex/DeepComplex/scripts/database_generating_scripts/hetero_std/test"
#local="/data/farhan/SoftwareTools/DeepComplex/DeepComplex/scripts/database_generating_scripts/hetero_std/all_dimer_atom_files"
#local="./test"
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
