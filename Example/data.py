import os, sys
os.environ['PYTHONHASHSEED']=str(2021)


import random
random.seed(2021) 

import numpy as np
np.random.seed(2021)


import MDAnalysis as mda
from CodeEntropy.FunctionCollection import EntropyFunctions as EF
from CodeEntropy.ClassCollection import DataContainer as DC

data_dir = os.path.dirname(os.path.abspath(__file__))
tprfile = os.path.join(data_dir,"data/md_A4_dna.tpr")
trrfile = os.path.join(data_dir,"data/md_A4_dna_xf.trr")
tScale = 1.0
fScale = 1.0
temper = 300.0 #K
u = mda.Universe(tprfile, trrfile)
thread = 8
axis_list = ["C5'", "C4'", "C3'"]
dataContainer = DC.DataContainer(u)
UA_entropyFF, UA_entropyTT, res_df = EF.compute_entropy_UA_level(
    arg_hostDataContainer = dataContainer,
    arg_outFile = None,
    arg_selector = 'all', 
    arg_moutFile = None,
    arg_nmdFile = None,
    arg_fScale = fScale,
    arg_tScale = tScale,
    arg_temper = temper,
    arg_verbose = 1,
    arg_axis_list = axis_list,
    arg_csv_out= "Atom_level_res_data.csv",
)
print (res_df)