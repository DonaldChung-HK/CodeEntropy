import os, sys
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import datetime as dt
from CodeEntropy.ClassCollection.PoseidonClass import Poseidon, Poseidon_mp

import pandas as pd
import dill as pickle
start = dt.datetime.now()
wd = os.path.dirname(os.path.abspath(__file__))
topo_file = os.path.join(wd,"data/1AKI_prod.tpr")
traj_file = os.path.join(wd,"data/1AKI_prod.trr")
u = mda.Universe(topo_file, traj_file)
# poseidon_object = Poseidon(container=u, start=2, end=12)
poseidon_object = Poseidon_mp(container=u, start=2, end=12, thread=2)
parsing = dt.datetime.now()
# print(poseidon_object.allMoleculeList)
result = poseidon_object.run_analysis(level_list = ['moleculeLevel'], verbose=False)
end = dt.datetime.now()
print (f"parsing taken = {parsing - start}")
print (f"time taken = {end - start}")