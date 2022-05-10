import os, sys
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader

from CodeEntropy.ClassCollection.PoseidonClass import Poseidon


def main():
    #hard coded for now
    wd = os.path.dirname(os.path.abspath(__file__))
    # loading files
    topo_file = os.path.join(wd,"data/1AKI_prod.tpr")
    traj_file = os.path.join(wd,"data/1AKI_prod.trr")
    # loading data into individual universe
    main = mda.Universe(topo_file, traj_file)
    poseidon_object = Poseidon(container=main, start=1, end=20, water=('SOL',), excludedResnames=("CL",))
    solventData, soluteData = poseidon_object.run_analysis(verbose=True)
    print("------------------------------SolventData------------------------------")
    print(solventData)
    print("------------------------------SoluteData------------------------------")
    print(soluteData)

if __name__ == '__main__':
    main()
