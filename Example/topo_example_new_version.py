# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 14:13:07 2022

@author: bmm66251
"""
import os, sys
import MDAnalysis as mda
from CodeEntropy.FunctionCollection import Utils
from CodeEntropy.FunctionCollection import EntropyFunctions
from CodeEntropy.IO import InputParser, Writer 
from CodeEntropy.ClassCollection import DataContainer as DC
from datetime import datetime

if __name__ == "__main__":
    startTime = datetime.now()
    ############## REPLACE INPUTS ##############
    wd = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(wd,"data/1AKI_ws.tpr")
    trrfile = os.path.join(wd,"data/1AKI_ws.trr")
    outfile = os.path.join(wd,"dna_mcc.out")
    outfile0_SC = os.path.join(wd,"dna_mcc_0_SC.out")
    outfile0_BB = os.path.join(wd,"dna_mcc_0_BB.out")
    outfile1_SC = os.path.join(wd,"dna_mcc_1_SC.out")
    outfile1_BB = os.path.join(wd,"dna_mcc_1_BB.out")
    tScale = 1
    fScale = 1
    temper = 300 #K
    
    u = mda.Universe(tprfile, trrfile)
    dataContainer = DC.DataContainer(u)
    
    result_entropy0_SC = EntropyFunctions.compute_topographical_entropy0_SC(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = outfile0_SC, 
        arg_verbose = 5
    )

    print(f"result_entropy0_SC = {result_entropy0_SC}")

    result_entropy0_BB = EntropyFunctions.compute_topographical_entropy0_BB(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = outfile0_BB, 
        arg_verbose = 5
    ) 

    print(f"result_entropy0_BB = {result_entropy0_BB}")


    result_entropy1_SC = EntropyFunctions.compute_topographical_entropy1_SC(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = outfile1_SC, 
        arg_verbose = 5
    )

    print(f"result_entropy1_SC= {result_entropy1_SC}")


    result_entropy1_BB = EntropyFunctions.compute_topographical_entropy1_BB(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = outfile1_BB, 
        arg_verbose = 5
    ) 
    print(f"result_entropy1_BB= {result_entropy1_BB}")
    
    result_entropy4 = EntropyFunctions.compute_topographical_entropy_method4(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = outfile1_SC, 
        arg_verbose = 5
    )

    print(f"result_entropy4 = {result_entropy4}")


    
    result_entropyAEM = EntropyFunctions.compute_topographical_entropy_AEM(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = outfile1_SC, 
        arg_verbose = 5
    )

    print(f"result_entropyAEM = {result_entropyAEM}")

    # Demanding computation if large amount of Dihedral
    result_entropy3 = EntropyFunctions.compute_topographical_entropy_method3(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = outfile1_BB, 
        arg_verbose = 5
    ) 
    print(f"result_entropy3 = {result_entropy3}")

    print(datetime.now() - startTime)
