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
import pandas as pd
import numpy as np
from datetime import datetime



if __name__ == "__main__":
    ############## REPLACE INPUTS ##############
    startTime = datetime.now()
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

    wm_entropyFF, wm_entropyTT = EF.compute_entropy_whole_molecule_level(
        arg_hostDataContainer = dataContainer,
        arg_outFile = None,
        arg_selector = "all", 
        arg_moutFile = None,
        arg_nmdFile = None,
        arg_fScale = fScale,
        arg_tScale = tScale,
        arg_temper = temper,
        arg_verbose = 5
    )

    print(f"wm_entropyFF = {wm_entropyFF}")
    print(f"wm_entropyTT = {wm_entropyTT}")

    res_entropyFF, res_entropyTT = EF.compute_entropy_residue_level(
        arg_hostDataContainer = dataContainer,
        arg_outFile = None,
        arg_selector = 'all', 
        arg_moutFile = None,
        arg_nmdFile = None,
        arg_fScale = fScale,
        arg_tScale = tScale,
        arg_temper = temper,
        arg_verbose = 5,
        arg_axis_list = axis_list,
    )
    print(res_entropyFF)
    print(res_entropyTT)


    UA_entropyFF, UA_entropyTT = EF.compute_entropy_UA_level(
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
        arg_csv_out= None,
    )
    print(f"UA_entropyFF = {UA_entropyFF}")
    print(f"UA_entropyTT = {UA_entropyTT}")

    UA_entropyFF, UA_entropyTT = EF.compute_entropy_UA_level_multiprocess(
        arg_hostDataContainer = dataContainer,
        arg_outFile = None,
        arg_selector = 'all', 
        arg_moutFile = None,
        arg_nmdFile = None,
        arg_fScale = fScale,
        arg_tScale = tScale,
        arg_temper = temper,
        arg_verbose = 1,
        arg_csv_out= None,
        arg_axis_list = axis_list,
        arg_thread= thread,
    )
    print(f"UA_entropyFF = {UA_entropyFF}")
    print(f"UA_entropyTT = {UA_entropyTT}")

    result_entropy0_SC = EF.compute_topographical_entropy0_SC(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )

    print(f"result_entropy0_SC = {result_entropy0_SC}")

    result_entropy0_BB = EF.compute_topographical_entropy0_BB(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    ) 

    print(f"result_entropy0_BB = {result_entropy0_BB}")

    

    result_entropy1_SC = EF.compute_topographical_entropy1_SC(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )

    print(f"result_entropy1_SC= {result_entropy1_SC}")


    result_entropy1_BB = EF.compute_topographical_entropy1_BB(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    ) 
    print(f"result_entropy1_BB= {result_entropy1_BB}")

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
        arg_outFile = None, 
        arg_verbose = 5
    )



    print(f"result_entropyAEM = {result_entropyAEM}")   
    result_entropy3 = EF.compute_topographical_entropy_method3(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = "entropy3.log", 
        arg_verbose = 5
    ) 
    print(f"result_entropy3 = {result_entropy3}")