from modeller import *
import os, sys

combined_pir_file=os.path.abspath(sys.argv[1])
outputfolder=os.path.abspath(sys.argv[2])+"/"
if not os.path.isdir(outputfolder): os.makedirs(outputfolder)

#Convert database to to .pir format extension .hdf5
log.verbose()
env = environ()
sdb = sequence_db(env)
sdb.convert(seq_database_file=combined_pir_file, seq_database_format='PIR',
            chains_list='ALL', minmax_db_seq_len=(30, 4000),
            clean_sequences=False, outfile=outputfolder+os.path.basename(combined_pir_file).split(".")[0]+".hdf5")