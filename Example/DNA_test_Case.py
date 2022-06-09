#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 10:30:58 2022

@author: donald
"""
import os, sys
import MDAnalysis as mda
from CodeEntropy.FunctionCollection import EntropyFunctions as EF
from CodeEntropy.ClassCollection import DataContainer as DC
from CodeEntropy.IO import MDAUniverseHelper as MDAHelper
import pandas as pd
import numpy as np
from datetime import datetime



if __name__ == "__main__":
    ############## REPLACE INPUTS ##############
    startTime = datetime.now()
    data_dir = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(data_dir,"data/md_A4_dna.tpr")
    trrfile = os.path.join(data_dir,"data/md_A4_dna_xf.trr")
    outfile = None
    tScale = 1.0
    fScale = 1.0
    temper = 300.0 #K
    u = mda.Universe(tprfile, trrfile)
    topo_outfile = None
    selection_string = None
    start = 0
    end = 101
    step = 2
    thread = 8
    axis_list = ["C5'", "C4'", "C3'"]
    if start == None:
        start = 0
    if end == None:
        end = len(u.trajectory)
    results_df = pd.DataFrame(columns=['Method and Level','Type', 'result'])
    dataContainer = DC.DataContainer(u)

    wm_entropyFF, wm_entropyTT = EF.compute_entropy_whole_molecule_level(
        arg_hostDataContainer = dataContainer,
        arg_outFile = None,
        arg_selector = "all", 
        arg_moutFile =  None,
        arg_nmdFile = None,
        arg_fScale = fScale,
        arg_tScale = tScale,
        arg_temper = None,
        arg_verbose = 0
    )
    
    res_entropyFF, res_entropyTT = EF.compute_entropy_residue_level(
        arg_hostDataContainer = dataContainer,
        arg_outFile = outfile,
        arg_selector = 'all', 
        arg_moutFile = 'ResidueLevel_matrix.out',
        arg_nmdFile = 'ResidueLevel_mode_spectra.out',
        arg_fScale = fScale,
        arg_tScale = tScale,
        arg_temper = temper,
        arg_verbose = 5,
        arg_axis_list = axis_list,
    )


    UA_entropyFF, UA_entropyTT = EF.compute_entropy_UA_level(
        arg_hostDataContainer = dataContainer,
        arg_outFile = outfile,
        arg_selector = 'all', 
        arg_moutFile = 'AtomLevel_matrix.out',
        arg_nmdFile = 'AtomLevel_mode_spectra.out',
        arg_fScale = fScale,
        arg_tScale = tScale,
        arg_temper = temper,
        arg_verbose = 1,
        arg_axis_list = axis_list,
        arg_csv_out= 'AtomLevel_bead_entropy.csv',
    )


    UA_entropyFF, UA_entropyTT = EF.compute_entropy_UA_level_multiprocess(
        arg_hostDataContainer = dataContainer,
        arg_outFile = outfile,
        arg_selector = 'all', 
        arg_moutFile = 'AtomLevel_matrix.out',
        arg_nmdFile = 'AtomLevel_mode_spectra.out',
        arg_fScale = fScale,
        arg_tScale = tScale,
        arg_temper = temper,
        arg_verbose = 1,
        arg_csv_out= 'AtomLevel_bead_entropy.csv',
        arg_axis_list = axis_list,
        arg_thread= thread,
    )


    result_entropy0_SC = EF.compute_topographical_entropy0_SC(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = topo_outfile, 
        arg_verbose = 5
    )


    result_entropy0_BB = EF.compute_topographical_entropy0_BB(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = topo_outfile, 
        arg_verbose = 5
    ) 

    

    result_entropy1_SC = EF.compute_topographical_entropy1_SC(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = topo_outfile, 
        arg_verbose = 5
    )


    result_entropy1_BB = EF.compute_topographical_entropy1_BB(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = topo_outfile, 
        arg_verbose = 5
    ) 


    # #work in progress
    # result_entropy4 = EntropyFunctions.compute_topographical_entropy_method4(
    #     arg_hostDataContainer = dataContainer, 
    #     arg_selector = "all",
    #     arg_outFile = outfile1_SC, 
    #     arg_verbose = 5
    # )

    # print(f"result_entropy4 = {result_entropy4}")


    
    result_entropyAEM = EF.compute_topographical_entropy_AEM(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = topo_outfile, 
        arg_verbose = 5
    )