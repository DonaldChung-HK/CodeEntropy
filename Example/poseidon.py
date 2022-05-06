import os, sys
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader

from CodeEntropy.poseidon.extractData import run_data_extractor 
from CodeEntropy.poseidon.analysis import run_poseidon_analysis

def load_data():
    #hard coded for now
    wd = os.path.dirname(os.path.abspath(__file__))
    # loading files
    topo_file = os.path.join(wd,"data/molecules.prmtop")
    traj_file = os.path.join(wd,"data/Trajectory_npt_1.data.gz")
    # remember to edit the format so that the header is "id mass x y z" otherwise MDAnalysis won't load the data 
    force_file = os.path.join(wd,"data/Forces_npt_1.data")
    energy_file = os.path.join(wd,"data/Energy_npt_1.data")
    # loading data into individual universe
    main = mda.Universe(topo_file, traj_file, atom_style='id type x y z' ,format="LAMMPSDUMP")
    force = mda.Universe(topo_file, force_file, atom_style='id type x y z' ,format="LAMMPSDUMP")
    energy = mda.Universe(topo_file, energy_file, atom_style='id type x y z' ,format="LAMMPSDUMP")
    # selection for accessing values
    select_atom = main.select_atoms('all')
    select_atom_force = force.select_atoms('all')
    select_atom_energy = energy.select_atoms('all')
    u2 = mda.Merge(select_atom)
    #loading values from universe
    coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(), select_atom).run().results['timeseries']
    forces = AnalysisFromFunction(lambda ag: ag.positions.copy(), select_atom_force).run().results['timeseries']
    energy = AnalysisFromFunction(lambda ag: ag.positions.copy(), select_atom_energy).run().results['timeseries']
    dimensions = AnalysisFromFunction(lambda ag: ag.dimensions.copy(), select_atom).run().results['timeseries']
    #merging it into one the energy is loaded into velocity since it is eeasier to selecting the value by frame/atom and there isn't a way to load custom data into a universe 
    u2.load_new(coordinates, format=MemoryReader, forces=forces, velocities=energy, dimensions=dimensions)
    return u2

def main():
    u2 = load_data() # this should drop all the intermediate values loaded which might benefit memory
    allMoleculeList = run_data_extractor.start(container=u2, start=1, end=20)
    run_poseidon_analysis.start(allMoleculeList)

if __name__ == '__main__':
    main()





