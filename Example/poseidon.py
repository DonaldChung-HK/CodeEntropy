import os, sys
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader

from CodeEntropy.poseidon.extractData import run_data_extractor 
from CodeEntropy.poseidon.analysis import run_poseidon_analysis

wd = os.path.dirname(os.path.abspath(__file__))
topo_file = os.path.join(wd,"data/molecules.prmtop")
traj_file = os.path.join(wd,"data/Trajectory_npt_1.data.gz")
force_file = os.path.join(wd,"data/Forces_npt_1.data")
energy_file = os.path.join(wd,"data/Energy_npt_1.data")
main = mda.Universe(topo_file, traj_file, atom_style='id type x y z' ,format="LAMMPSDUMP")
force = mda.Universe(topo_file, force_file, atom_style='id type x y z' ,format="LAMMPSDUMP")
energy = mda.Universe(topo_file, energy_file, atom_style='id type x y z' ,format="LAMMPSDUMP")
select_atom = main.select_atoms('all')
select_atom_force = force.select_atoms('all')
select_atom_energy = energy.select_atoms('all')
u2 = mda.Merge(select_atom)
coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(), select_atom).run().results['timeseries']
forces = AnalysisFromFunction(lambda ag: ag.positions.copy(), select_atom_force).run().results['timeseries']
energy = AnalysisFromFunction(lambda ag: ag.positions.copy(), select_atom_energy).run().results['timeseries']
dimensions = AnalysisFromFunction(lambda ag: ag.dimensions.copy(), select_atom).run().results['timeseries']
u2.load_new(coordinates, format=MemoryReader, forces=forces, velocities=energy, dimensions=dimensions)

allMoleculeList = run_data_extractor.start(start=1, end=20, step=1, container=u2)
run_poseidon_analysis.start(allMoleculeList)