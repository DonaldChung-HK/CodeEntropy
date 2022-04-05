import os, sys
import MDAnalysis as mda
from CodeEntropy.FunctionCollection import EntropyFunctions as EF
from CodeEntropy.ClassCollection import DataContainer as DC

if __name__ == "__main__":
    ############## REPLACE INPUTS ##############
    wd = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(wd,"data/1AKI_ws.tpr")
    trrfile = os.path.join(wd,"data/1AKI_ws.trr")
    outfile = os.path.join(wd,"1AKI_mcc.out")
    tScale = 1
    fScale = 1
    temper = 300 #K
    u = mda.Universe(tprfile, trrfile)
    dataContainer = DC.DataContainer(u)

    wm_entropyFF, wm_entropyTT = EF.compute_entropy_whole_molecule_level(
        arg_hostDataContainer = dataContainer,
        arg_outFile = outfile,
        arg_selector = "all", 
        arg_moutFile = None,
        arg_nmdFile = None,
        arg_fScale = 1,
        arg_tScale = 1,
        arg_temper = 300,
        arg_verbose = 3
    )

    res_entropyFF, res_entropyTT = EF.compute_entropy_residue_level(
        arg_hostDataContainer = dataContainer,
        arg_outFile = outfile,
        arg_selector = "all", 
        arg_moutFile = None,
        arg_nmdFile = None,
        arg_fScale = 1,
        arg_tScale = 1,
        arg_temper = 300,
        arg_verbose = 3
    )
    
    UA_entropyFF, UA_entropyTT = EF.compute_entropy_UA_level(
        arg_hostDataContainer = dataContainer,
        arg_outFile = outfile,
        arg_selector = "all", 
        arg_moutFile = None,
        arg_nmdFile = None,
        arg_fScale = 1,
        arg_tScale = 1,
        arg_temper = 300,
        arg_verbose = 3
    )