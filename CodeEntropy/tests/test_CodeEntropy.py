import os, sys
import MDAnalysis as mda
from CodeEntropy.FunctionCollection import EntropyFunctions as EF
from CodeEntropy.ClassCollection import DataContainer as DC

import pytest

def test_CodeEntropy_whole_molecule():
    """test for computing entropy at whole molecule level"""
    data_dir = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(data_dir,"data/md_A4_dna.tpr")
    trrfile = os.path.join(data_dir,"data/md_A4_dna_xf.trr")
    tScale = 1.0
    fScale = 1.0
    temper = 300.0 #K
    u = mda.Universe(tprfile, trrfile)
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
    assert wm_entropyFF == pytest.approx(47.230063169278864)
    assert wm_entropyTT == pytest.approx(44.35619230374928)

def test_CodeEntropy_res_level():
    """test for computing entropy at residue level"""
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
    assert res_entropyFF == pytest.approx(194.93511803752938)
    assert res_entropyTT == pytest.approx(285.7213867860228)

def test_CodeEntropy_united_atom_level():
    """test for computing entropy at united atom level"""
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

    assert UA_entropyFF == pytest.approx(1439.501746003283)
    assert UA_entropyTT == pytest.approx(107.67610435497281)

def test_CodeEntropy_united_atom_level_multiprocess():
    """test for computing entropy at united atom level with multiprocess"""
    data_dir = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(data_dir,"data/md_A4_dna.tpr")
    trrfile = os.path.join(data_dir,"data/md_A4_dna_xf.trr")
    tScale = 1.0
    fScale = 1.0
    temper = 300.0 #K
    u = mda.Universe(tprfile, trrfile)
    thread = 2
    axis_list = ["C5'", "C4'", "C3'"]
    dataContainer = DC.DataContainer(u)
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
    assert UA_entropyFF == pytest.approx(1439.5017460032827)
    assert UA_entropyTT == pytest.approx(107.67610435497281)