import os, sys
import MDAnalysis as mda
import datetime as dt
from CodeEntropy.ClassCollection.PoseidonClass import Poseidon
import pandas as pd
import pytest
from numpy import array, testing

@pytest.fixture(scope='session')
def wd():
    wd = os.path.dirname(os.path.abspath(__file__))
    return wd
@pytest.fixture(scope='session')
def poseidon_object(wd):
    topo_file = os.path.join(wd,"data/poseidon_example.prmtop")
    traj_file = os.path.join(wd,"data/poseidon_example.trr")
    u = mda.Universe(topo_file, traj_file)
    poseidon_object = Poseidon(container=u, start=2, end=12)
    return poseidon_object

parser_ref_index = [
    (0, 9),
    (1, 916), 
    (2, 4584), 
    (3, 2093), 
    (4, 8641), 
    (5, 5221)
]
@pytest.fixture(scope='session')
def parser_ref():
    parser_ref = [
        [2,
        'CZ',
        21,
        'ARG',
        2,
        [3, 0],
        [6, 8, 10, 13, 16, 19, 21, 22, 25, 28, 29],
        ['O', 'WAT', 240, 3.307],
        1,
        (['0', 'X'], [], []),
        None,
        None,
        None,
        array([[ 14.508, -29.179, -13.636],
                [-29.179,  58.689,  27.426],
                [-13.636,  27.426,  12.816]]),
        array([[0., 0., 0.],
                [0., 0., 0.],
                [0., 0., 0.]]),
        None,
        None,
        None,
        [['CD', 2, 'ARG'], ['O', 240, 'WAT']]],
        [3,
        'CH3',
        1,
        'ACE',
        1,
        [1, 3],
        [1, 4, 5],
        ['N', 'ARG', 2, 2.428],
        1,
        (['0', '1', 'WAT_O', 'WAT_O', 'WAT_O', 'WAT_O', 'WAT_O'], [], []),
        array([[ 0.258, -0.   ,  0.206],
                [-0.   ,  0.   , -0.   ],
                [ 0.206, -0.   ,  0.165]]),
        array([[ 4.133,  1.435, -6.39 ],
                [ 1.435,  0.498, -2.219],
                [-6.39 , -2.219,  9.882]]),
        None,
        array([[ 2.000e-03,  6.400e-02, -1.420e-01],
                [ 6.400e-02,  1.845e+00, -4.090e+00],
                [-1.420e-01, -4.090e+00,  9.066e+00]]),
        array([[ 0.657,  5.304,  4.153],
                [ 5.304, 42.808, 33.517],
                [ 4.153, 33.517, 26.243]]),
        [],
        array([[  0.,   0.,  -1.,   0.,  -0.,   1.,  -0.,  -0.,  -0.],
                [  0.,   7., -16.,  12.,  -2.,  32.,  -6.,  -6.,  -8.],
                [ -1., -16.,  36., -27.,   4., -70.,  14.,  13.,  17.],
                [  0.,  12., -27.,  20.,  -3.,  51., -10.,  -9., -13.],
                [ -0.,  -2.,   4.,  -3.,   0.,  -7.,   1.,   1.,   2.],
                [  1.,  32., -70.,  51.,  -7., 134., -26., -24., -33.],
                [ -0.,  -6.,  14., -10.,   1., -26.,   5.,   5.,   6.],
                [ -0.,  -6.,  13.,  -9.,   1., -24.,   5.,   4.,   6.],
                [ -0.,  -8.,  17., -13.,   2., -33.,   6.,   6.,   8.]]),
        array([[ 0.657,  5.304,  4.153,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                0.   ],
                [ 5.304, 42.808, 33.517,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                0.   ],
                [ 4.153, 33.517, 26.243,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                0.   ],
                [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                0.   ],
                [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                0.   ],
                [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                0.   ],
                [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                0.   ],
                [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                0.   ],
                [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                0.   ]]),
        [['O', 1, 'ACE'],
        ['N', 2, 'ARG'],
        ['O', 817, 'WAT'],
        ['O', 61, 'WAT'],
        ['O', 539, 'WAT'],
        ['O', 341, 'WAT'],
        ['O', 315, 'WAT']]],
        [7,
        'CA',
        8,
        'ARG',
        2,
        [3, 1],
        [6, 8, 10, 13, 16, 19, 21, 22, 25, 28, 29],
        ['C', 'ACE', 1, 2.484],
        1,
        (['0', 'NME_N', 'WAT_O', 'X', 'X'], [], []),
        None,
        None,
        None,
        array([[ 5.057, 16.631, 18.57 ],
                [16.631, 54.691, 61.068],
                [18.57 , 61.068, 68.188]]),
        array([[17.433, 14.215, 33.262],
                [14.215, 11.59 , 27.121],
                [33.262, 27.121, 63.463]]),
        None,
        None,
        None,
        [['C', 1, 'ACE'],
        ['O', 2, 'ARG'],
        ['N', 3, 'NME'],
        ['CG', 2, 'ARG'],
        ['O', 422, 'WAT']]],
        [4,
        'O',
        771,
        'WAT',
        249,
        [0, 2],
        [771],
        [],
        2,
        ([], [], []),
        array([[ 0.009,  0.027, -0.004],
                [ 0.027,  0.076, -0.01 ],
                [-0.004, -0.01 ,  0.001]]),
        array([[ 0.459,  0.477, -0.352],
                [ 0.477,  0.496, -0.365],
                [-0.352, -0.365,  0.269]]),
        None,
        array([[ 0.009,  0.027, -0.004],
                [ 0.027,  0.076, -0.01 ],
                [-0.004, -0.01 ,  0.001]]),
        array([[ 0.459,  0.477, -0.352],
                [ 0.477,  0.496, -0.365],
                [-0.352, -0.365,  0.269]]),
        None,
        array([[ 0.009,  0.027, -0.004],
                [ 0.027,  0.076, -0.01 ],
                [-0.004, -0.01 ,  0.001]]),
        array([[ 0.459,  0.477, -0.352],
                [ 0.477,  0.496, -0.365],
                [-0.352, -0.365,  0.269]]),
        []],
        [11,
        'O',
        1179,
        'WAT',
        385,
        [0, 2],
        [1179],
        [],
        2,
        ([], [], []),
        array([[0.   , 0.   , 0.   ],
                [0.   , 0.351, 0.102],
                [0.   , 0.102, 0.03 ]]),
        array([[ 0.77 ,  0.908, -1.873],
                [ 0.908,  1.07 , -2.208],
                [-1.873, -2.208,  4.557]]),
        None,
        array([[0.   , 0.   , 0.   ],
                [0.   , 0.351, 0.102],
                [0.   , 0.102, 0.03 ]]),
        array([[ 0.77 ,  0.908, -1.873],
                [ 0.908,  1.07 , -2.208],
                [-1.873, -2.208,  4.557]]),
        None,
        array([[0.   , 0.   , 0.   ],
                [0.   , 0.351, 0.102],
                [0.   , 0.102, 0.03 ]]),
        array([[ 0.77 ,  0.908, -1.873],
                [ 0.908,  1.07 , -2.208],
                [-1.873, -2.208,  4.557]]),
        []],
        [7,
        'O',
        1911,
        'WAT',
        629,
        [0, 2],
        [1911],
        [],
        2,
        ([], [], []),
        array([[0.024, 0.05 , 0.223],
                [0.05 , 0.103, 0.459],
                [0.223, 0.459, 2.037]]),
        array([[ 2.0345e+01, -4.8700e-01, -1.1815e+01],
                [-4.8700e-01,  1.2000e-02,  2.8300e-01],
                [-1.1815e+01,  2.8300e-01,  6.8610e+00]]),
        None,
        array([[0.024, 0.05 , 0.223],
                [0.05 , 0.103, 0.459],
                [0.223, 0.459, 2.037]]),
        array([[ 2.0345e+01, -4.8700e-01, -1.1815e+01],
                [-4.8700e-01,  1.2000e-02,  2.8300e-01],
                [-1.1815e+01,  2.8300e-01,  6.8610e+00]]),
        None,
        array([[0.024, 0.05 , 0.223],
                [0.05 , 0.103, 0.459],
                [0.223, 0.459, 2.037]]),
        array([[ 2.0345e+01, -4.8700e-01, -1.1815e+01],
                [-4.8700e-01,  1.2000e-02,  2.8300e-01],
                [-1.1815e+01,  2.8300e-01,  6.8610e+00]]),
        []]
    ]
    return parser_ref

# the reason for testing individually is to allow identification of which value went wrong
@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_frame(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 0 which store the frame of the atom info
    """
    assert poseidon_object.allMoleculeList[test_id][0] == parser_ref[ref_index][0]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_atom_name(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 1 which store the atom name
    """
    assert poseidon_object.allMoleculeList[test_id][1] == parser_ref[ref_index][1]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_atom_num(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 2 which store the atom number
    """
    assert poseidon_object.allMoleculeList[test_id][2] == parser_ref[ref_index][2]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_resname(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 3 which store the resname
    """
    assert poseidon_object.allMoleculeList[test_id][3] == parser_ref[ref_index][3]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_resid(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 4 which store the resid
    """
    assert poseidon_object.allMoleculeList[test_id][4] == parser_ref[ref_index][4]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_bondedUA_H(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 5 which store the bondedUA_H which is a tuple of [num_bondedUAs, num_bonded Hs]
    """
    assert poseidon_object.allMoleculeList[test_id][5] == parser_ref[ref_index][5]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_molecule_atomNums(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 6 which store the list of bonded UA atom nums
    """
    assert poseidon_object.allMoleculeList[test_id][6] == parser_ref[ref_index][6]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_nearestInfo(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 7 which store the list of nearest non like neighbour in a tuple of (atom_name, resname, resid, distance)
    """
    assert poseidon_object.allMoleculeList[test_id][7] == parser_ref[ref_index][7]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_hydrationShellRAD_solute(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 8 which store the list of nearest solute molecule centric
    """
    assert poseidon_object.allMoleculeList[test_id][8] == parser_ref[ref_index][8]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_RAD_shell_ranked(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 9 which store the list of 3 tuple (RAD shell ranked, acceptors ranked, donors ranked by RAD)
    """
    assert poseidon_object.allMoleculeList[test_id][9] == parser_ref[ref_index][9]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_MweightedForces(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 10 which store the MweightedForces
    """
    #numpy can't compare None type object
    if not parser_ref[ref_index][10] is None:
        testing.assert_array_almost_equal(poseidon_object.allMoleculeList[test_id][10], parser_ref[ref_index][10])
    else:
        assert poseidon_object.allMoleculeList[test_id][10] == parser_ref[ref_index][10]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_MweightedTorques(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 11 which store the MweightedTorques
    """
    #numpy can't compare None type object
    if not parser_ref[ref_index][11] is None:
        testing.assert_array_almost_equal(poseidon_object.allMoleculeList[test_id][11], parser_ref[ref_index][11])
    else:
        assert poseidon_object.allMoleculeList[test_id][11] == parser_ref[ref_index][11]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_UA_PKenergy(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 12 which store the UA_PKenergy should all be None as energy is out of scope
    """
    assert poseidon_object.allMoleculeList[test_id][12] == parser_ref[ref_index][12]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_UAweightedForces(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 13 which store the UAweightedForces
    """
    #numpy can't compare None type object
    if not parser_ref[ref_index][13] is None:
        testing.assert_array_almost_equal(poseidon_object.allMoleculeList[test_id][13], parser_ref[ref_index][13])
    else:
        assert poseidon_object.allMoleculeList[test_id][13] == parser_ref[ref_index][13]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_UAweightedTorques(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 14 which store the UAweightedTorques
    """
    #numpy can't compare None type object
    if not parser_ref[ref_index][14] is None:
        testing.assert_array_almost_equal(poseidon_object.allMoleculeList[test_id][14], parser_ref[ref_index][14])
    else:
        assert poseidon_object.allMoleculeList[test_id][14] == parser_ref[ref_index][14]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_dihedral_phi_type(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 15 which store the dihedral_phi_type a list of unique dihedral angles, in degrees
    """
    assert poseidon_object.allMoleculeList[test_id][15] == parser_ref[ref_index][15]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_molecule_UA_Fs(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 16 which store the molecule_UA_Fs
    """
    #numpy can't compare None type object
    if not parser_ref[ref_index][16] is None:
        testing.assert_array_almost_equal(poseidon_object.allMoleculeList[test_id][16], parser_ref[ref_index][16])
    else:
        assert poseidon_object.allMoleculeList[test_id][16] == parser_ref[ref_index][16]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_molecule_UA_Ts(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 17 which store the molecule_UA_Ts
    """
    #numpy can't compare None type object
    if not parser_ref[ref_index][17] is None:
        testing.assert_array_almost_equal(poseidon_object.allMoleculeList[test_id][17], parser_ref[ref_index][17])
    else:
        assert poseidon_object.allMoleculeList[test_id][17] == parser_ref[ref_index][17]

@pytest.mark.parametrize("ref_index,test_id",
    parser_ref_index
)
def test_poseidon_parser_RAD_nAtoms(poseidon_object, parser_ref, ref_index,test_id):
    """
    testing the return of allMoleculeList index 18 which store the RAD_nAtoms
    """
    assert poseidon_object.allMoleculeList[test_id][18] == parser_ref[ref_index][18]




def test_poseidon_moleculeLevel(poseidon_object, wd):
    """
    testing result of poseidon at moleculeLevel
    """
    result = poseidon_object.run_analysis(level_list = ['moleculeLevel'], verbose=False)
    data_dir_solute = os.path.join(wd,"data/soluteVariables10.0EE_moleculeLevel.csv")
    solute_ref = pd.read_csv(data_dir_solute, na_values="nan")
    pd.testing.assert_frame_equal(solute_ref, result["moleculeLevel"]["soluteData"], check_dtype=False)
    data_dir_solvent = os.path.join(wd,"data/solventVariables10.0EE_moleculeLevel.csv")
    solvent_ref = pd.read_csv(data_dir_solvent, na_values="nan")
    pd.testing.assert_frame_equal(solvent_ref, result["moleculeLevel"]["solventData"], check_dtype=False)
    
def test_poseidon_residLevel_resname(poseidon_object, wd):
    """
    testing result of poseidon at residLevel_resname level
    """
    result = poseidon_object.run_analysis(level_list = ['residLevel_resname'], verbose=False)
    data_dir_solute = os.path.join(wd,"data/soluteVariables10.0EE_residLevel_resname.csv")
    solute_ref = pd.read_csv(data_dir_solute, na_values="nan")
    pd.testing.assert_frame_equal(solute_ref, result["residLevel_resname"]["soluteData"], check_dtype=False)
    data_dir_solvent = os.path.join(wd,"data/solventVariables10.0EE_residLevel_resname.csv")
    solvent_ref = pd.read_csv(data_dir_solvent, na_values="nan")
    pd.testing.assert_frame_equal(solvent_ref, result["residLevel_resname"]["solventData"], check_dtype=False)
    data_dir_contact = os.path.join(wd,"data/resid_contact_matrix_residLevel_resname.csv")
    contact_ref = pd.read_csv(data_dir_contact, na_values="nan")
    pd.testing.assert_frame_equal(contact_ref, result["residLevel_resname"]["contactMatrix"], check_dtype=False)

def test_poseidon_atomLevel(poseidon_object, wd):
    """
    testing result of poseidon at atomLevel
    """
    result = poseidon_object.run_analysis(level_list = ['atomLevel'], verbose=False)
    data_dir_solute = os.path.join(wd,"data/soluteVariables10.0EE_atomLevel.csv")
    solute_ref = pd.read_csv(data_dir_solute, na_values="nan")
    pd.testing.assert_frame_equal(solute_ref, result["atomLevel"]["soluteData"], check_dtype=False)
    data_dir_solvent = os.path.join(wd,"data/solventVariables10.0EE_atomLevel.csv")
    solvent_ref = pd.read_csv(data_dir_solvent, na_values="nan")
    pd.testing.assert_frame_equal(solvent_ref, result["atomLevel"]["solventData"], check_dtype=False)
    data_dir_contact = os.path.join(wd,"data/resid_contact_matrix_atomLevel.csv")
    contact_ref = pd.read_csv(data_dir_contact, na_values="nan")
    pd.testing.assert_frame_equal(contact_ref, result["atomLevel"]["contactMatrix"], check_dtype=False)

def test_poseidon_soluteContacts(poseidon_object, wd):
    """
    testing result of poseidon at soluteContacts Level
    """
    result = poseidon_object.run_analysis(level_list = ['soluteContacts'], verbose=False)
    data_dir_solute = os.path.join(wd,"data/soluteVariables10.0EE_soluteContacts.csv")
    solute_ref = pd.read_csv(data_dir_solute, na_values="nan")
    pd.testing.assert_frame_equal(solute_ref, result["soluteContacts"]["soluteData"], check_dtype=False)
    data_dir_solvent = os.path.join(wd,"data/solventVariables10.0EE_soluteContacts.csv")
    solvent_ref = pd.read_csv(data_dir_solvent, na_values="nan")
    pd.testing.assert_frame_equal(solvent_ref, result["soluteContacts"]["solventData"], check_dtype=False)