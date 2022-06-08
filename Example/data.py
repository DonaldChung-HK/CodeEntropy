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
tprfile = os.path.join(data_dir,"data/1AKI_ws.tpr")
trrfile = os.path.join(data_dir,"data/1AKI_ws.trr")
u = mda.Universe(tprfile, trrfile)
dataContainer = DC.DataContainer(u)
result_entropy0_SC = EF.compute_topographical_entropy_AEM(
    arg_hostDataContainer = dataContainer, 
    arg_selector = "all",
    arg_outFile = None, 
    arg_verbose = 5
)

print(f"result_entropy0_SC = {result_entropy0_SC}")