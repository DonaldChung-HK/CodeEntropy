##########################################################################################
# !!!!!!!!!!!The 1AKI_prod.trr trajectory is too big and can't be included in the repo
##########################################################################################

import os, sys
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import pandas as pd
from CodeEntropy.ClassCollection.PoseidonClass import Poseidon, Poseidon_mp


def main():
    #hard coded for now
    wd = os.path.dirname(os.path.abspath(__file__))
    # loading files
    topo_file = os.path.join(wd,"data/1AKI_prod.tpr")
    traj_file = os.path.join(wd,"data/1AKI_prod.trr")
    # loading data into individual universe
    main = mda.Universe(topo_file, traj_file)
    poseidon_object = Poseidon_mp(container=main, start=1, end=30, step=2,water=('SOL',), excludedResnames=("CL",), verbose=True)
    poseidon_object.run_analysis(level_list=['moleculeLevel', 'residLevel_resname', 'atomLevel', 'soluteContact'], verbose=True)

if __name__ == '__main__':
    main()
