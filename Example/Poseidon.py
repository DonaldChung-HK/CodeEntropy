import os, sys
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import datetime as dt
from CodeEntropy.ClassCollection.PoseidonClass import Poseidon, Poseidon_mp

import pandas as pd

wd = os.path.dirname(os.path.abspath(__file__))
topo_file = os.path.join(wd,"data/poseidon_example.prmtop")
traj_file = os.path.join(wd,"data/poseidon_example.trr")
u = mda.Universe(topo_file, traj_file)
poseidon_object = Poseidon_mp(container=u, start=2, end=12)
result = poseidon_object.run_analysis(level_list = ['moleculeLevel'], verbose=False)