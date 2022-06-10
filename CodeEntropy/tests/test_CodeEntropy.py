import os, sys
import MDAnalysis as mda
from CodeEntropy.FunctionCollection import EntropyFunctions as EF
from CodeEntropy.ClassCollection import DataContainer as DC
import numpy as np
import pandas as pd
import pytest


@pytest.fixture(scope='session')
def data_dir():
    data_dir = os.path.dirname(os.path.abspath(__file__))
    return data_dir
@pytest.fixture(scope='session')
def dataContainer_DNA(data_dir):
    tprfile = os.path.join(data_dir,"data/md_A4_dna.tpr")
    trrfile = os.path.join(data_dir,"data/md_A4_dna_xf.trr")
    u = mda.Universe(tprfile, trrfile)
    dataContainer_DNA = DC.DataContainer(u)
    return dataContainer_DNA

@pytest.fixture(scope='session')
def dataContainer_protein(data_dir):
    tprfile = os.path.join(data_dir,"data/1AKI_ws.tpr")
    trrfile = os.path.join(data_dir,"data/1AKI_ws.trr")
    u = mda.Universe(tprfile, trrfile)
    dataContainer_protein = DC.DataContainer(u)
    return dataContainer_protein

@pytest.fixture(scope='session')
def tScale():
    tScale = 1.0
    return tScale
@pytest.fixture(scope='session')
def fScale():
    fScale = 1.0
    return fScale
@pytest.fixture(scope='session')
def temper():
    temper = 300.0 #K
    return temper
@pytest.fixture(scope='session')
def axis_list():
    axis_list = ["C5'", "C4'", "C3'"]
    return axis_list
@pytest.fixture(scope='session')
def thread():
    tScale = 2
    return 2


def test_CodeEntropy_whole_molecule(dataContainer_DNA, tScale, fScale, temper):
    """test for computing entropy at whole molecule level"""
    wm_entropyFF, wm_entropyTT = EF.compute_entropy_whole_molecule_level(
        arg_hostDataContainer = dataContainer_DNA,
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

def test_CodeEntropy_res_level(dataContainer_DNA, tScale, fScale, temper, axis_list):
    """test for computing entropy at residue level"""
    res_entropyFF, res_entropyTT = EF.compute_entropy_residue_level(
        arg_hostDataContainer = dataContainer_DNA,
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


def test_CodeEntropy_united_atom_level(dataContainer_DNA, tScale, fScale, temper, axis_list, data_dir):
    """test for computing entropy at united atom level"""
    UA_entropyFF, UA_entropyTT, res_df = EF.compute_entropy_UA_level(
        arg_hostDataContainer = dataContainer_DNA,
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
    csv_dir = os.path.join(data_dir,"data/Atom_level_res_data.csv")
    res_ref = pd.read_csv(csv_dir, na_values="nan")
    pd.testing.assert_frame_equal(res_ref, res_df, check_dtype=False)

def test_CodeEntropy_united_atom_level_multiprocess(dataContainer_DNA, tScale, fScale, temper, axis_list, thread, data_dir):
    """test for computing entropy at united atom level with multiprocess"""
    UA_entropyFF, UA_entropyTT, res_df = EF.compute_entropy_UA_level_multiprocess(
        arg_hostDataContainer = dataContainer_DNA,
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
    csv_dir = os.path.join(data_dir,"data/Atom_level_res_data.csv")
    res_ref = pd.read_csv(csv_dir, na_values="nan")
    pd.testing.assert_frame_equal(res_ref, res_df, check_dtype=False)

def test_CodeEntropy_topo_method5(dataContainer_protein):
    """test for computing topographical entroy using method 5 AEM method. AEM is random so it is impossible to make test accurately"""
    result_entropyAEM = EF.compute_topographical_entropy_AEM(
        arg_hostDataContainer = dataContainer_protein, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )
    assert 220.0 < result_entropyAEM < 230.0

def test_CodeEntropy_topo_method3(dataContainer_DNA):
    """test for computing topographical entroy using method 3 Correlation density method"""
    result_entropy3 = EF.compute_topographical_entropy_method3(
        arg_hostDataContainer = dataContainer_DNA, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    ) 
    result = np.sum(result_entropy3)
    reference = 211.99999999999983
    assert reference == pytest.approx(result)

def test_CodeEntropy_topo_method0_SC(dataContainer_protein):
    """test for computing topographical entroy using method 3 Correlation density method"""
    result_entropy0_SC = EF.compute_topographical_entropy0_SC(
        arg_hostDataContainer = dataContainer_protein, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )

    reference = 621.2072742842519
    assert reference == pytest.approx(result_entropy0_SC)

def test_CodeEntropy_topo_method0_BB(dataContainer_protein):
    """test for computing topographical entroy using method 3 Correlation density method"""
    result_entropy0_BB = EF.compute_topographical_entropy0_BB(
        arg_hostDataContainer = dataContainer_protein, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )

    reference = 728.8710286898643
    assert reference == pytest.approx(result_entropy0_BB)

def test_CodeEntropy_topo_method1_SC(dataContainer_protein):
    """test for computing topographical entroy using method 1 some degree of randomness so a range is used"""
    result_entropy1_SC = EF.compute_topographical_entropy1_SC(
        arg_hostDataContainer = dataContainer_protein, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )

    reference = 390.587
    assert 388.0 <= result_entropy1_SC <= 395.0

def test_CodeEntropy_topo_method1_BB(dataContainer_protein):
    """test for computing topographical entroy using method 1"""
    result_entropy1_BB = EF.compute_topographical_entropy1_BB(
        arg_hostDataContainer = dataContainer_protein, 
        arg_selector = "all",
        arg_outFile = None, 
        arg_verbose = 5
    )

    reference = 38.37223423319515
    assert reference == pytest.approx(result_entropy1_BB)
