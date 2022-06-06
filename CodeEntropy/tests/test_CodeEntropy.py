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
