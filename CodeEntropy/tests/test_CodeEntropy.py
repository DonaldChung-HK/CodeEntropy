import os, sys
import MDAnalysis as mda
from CodeEntropy.FunctionCollection import EntropyFunctions as EF
from CodeEntropy.ClassCollection import DataContainer as DC
import numpy as np
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

def test_CodeEntropy_topo_method5():
    """test for computing topographical entroy using method 5 AEM method. AEM is random so it is impossible to make test accurately"""
    data_dir = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(data_dir,"data/1AKI_ws.tpr")
    trrfile = os.path.join(data_dir,"data/1AKI_ws.trr")
    u = mda.Universe(tprfile, trrfile)
    dataContainer = DC.DataContainer(u)
    result_entropyAEM = EF.compute_topographical_entropy_AEM(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )
    assert 220.0 < result_entropyAEM < 230.0

def test_CodeEntropy_topo_method3():
    """test for computing topographical entroy using method 3 Correlation density method"""
    data_dir = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(data_dir,"data/md_A4_dna.tpr")
    trrfile = os.path.join(data_dir,"data/md_A4_dna_xf.trr")
    u = mda.Universe(tprfile, trrfile)
    dataContainer = DC.DataContainer(u)
    result_entropy3 = EF.compute_topographical_entropy_method3(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    ) 
    result = np.sum(result_entropy3)
    reference = 211.99999999999983
    assert reference == pytest.approx(result)

def test_CodeEntropy_topo_method0_SC():
    """test for computing topographical entroy using method 3 Correlation density method"""
    data_dir = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(data_dir,"data/1AKI_ws.tpr")
    trrfile = os.path.join(data_dir,"data/1AKI_ws.trr")
    u = mda.Universe(tprfile, trrfile)
    dataContainer = DC.DataContainer(u)
    result_entropy0_SC = EF.compute_topographical_entropy0_SC(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )

    reference = 621.2072742842519
    assert reference == pytest.approx(result_entropy0_SC)

def test_CodeEntropy_topo_method0_BB():
    """test for computing topographical entroy using method 3 Correlation density method"""
    data_dir = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(data_dir,"data/1AKI_ws.tpr")
    trrfile = os.path.join(data_dir,"data/1AKI_ws.trr")
    u = mda.Universe(tprfile, trrfile)
    dataContainer = DC.DataContainer(u)
    result_entropy0_BB = EF.compute_topographical_entropy0_BB(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )

    reference = 728.8710286898643
    assert reference == pytest.approx(result_entropy0_BB)

def test_CodeEntropy_topo_method1_SC():
    """test for computing topographical entroy using method 1 some degree of randomness so a range is used"""
    data_dir = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(data_dir,"data/1AKI_ws.tpr")
    trrfile = os.path.join(data_dir,"data/1AKI_ws.trr")
    u = mda.Universe(tprfile, trrfile)
    dataContainer = DC.DataContainer(u)
    result_entropy1_SC = EF.compute_topographical_entropy1_SC(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )

    reference = 390.587
    assert 390.0 <= result_entropy1_SC <= 391

def test_CodeEntropy_topo_method1_BB():
    """test for computing topographical entroy using method 1"""
    data_dir = os.path.dirname(os.path.abspath(__file__))
    tprfile = os.path.join(data_dir,"data/1AKI_ws.tpr")
    trrfile = os.path.join(data_dir,"data/1AKI_ws.trr")
    u = mda.Universe(tprfile, trrfile)
    dataContainer = DC.DataContainer(u)
    result_entropy1_BB = EF.compute_topographical_entropy1_BB(
        arg_hostDataContainer = dataContainer, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )

    reference = 38.37223423319515
    assert reference == pytest.approx(result_entropy1_BB)
