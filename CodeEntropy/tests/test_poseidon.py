import os, sys
import MDAnalysis as mda
import datetime as dt
from CodeEntropy.ClassCollection.PoseidonClass import Poseidon
import pandas as pd

def test_poseidon_moleculeLevel():
    start = dt.datetime.now()
    wd = os.path.dirname(os.path.abspath(__file__))
    topo_file = os.path.join(wd,"data/poseidon_example.prmtop")
    traj_file = os.path.join(wd,"data/poseidon_example.trr")
    u = mda.Universe(topo_file, traj_file)
    poseidon_object = Poseidon(container=u, start=2, end=12)
    result = poseidon_object.run_analysis(level_list = ['moleculeLevel'], verbose=False)
    data_dir_solute = os.path.join(wd,"data/soluteVariables10.0EE_moleculeLevel.csv")
    solute_ref = pd.read_csv(data_dir_solute, na_values="nan")
    pd.testing.assert_frame_equal(solute_ref, result["moleculeLevel"]["soluteData"], check_dtype=False)
    data_dir_solvent = os.path.join(wd,"data/solventVariables10.0EE_moleculeLevel.csv")
    solvent_ref = pd.read_csv(data_dir_solvent, na_values="nan")
    pd.testing.assert_frame_equal(solvent_ref, result["moleculeLevel"]["solventData"], check_dtype=False)
    
