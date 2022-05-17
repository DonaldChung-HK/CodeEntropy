import os, sys
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import datetime as dt
from CodeEntropy.ClassCollection.PoseidonClass import Poseidon

def load_data():
    #hard coded for now
    wd = os.path.dirname(os.path.abspath(__file__))
    # loading files
    topo_file = os.path.join(wd,"data/molecules.prmtop")
    traj_file = os.path.join(wd,"data/Trajectory_npt_1.data.gz")
    # remember to edit the format so that the header is "id mass x y z" otherwise MDAnalysis won't load the data 
    force_file = os.path.join(wd,"data/Forces_npt_1.data")
    # loading data into individual universe
    main = mda.Universe(topo_file, traj_file, atom_style='id type x y z' ,format="LAMMPSDUMP")
    force = mda.Universe(topo_file, force_file, atom_style='id type x y z' ,format="LAMMPSDUMP")
    # selection for accessing values
    select_atom = main.select_atoms('all')
    select_atom_force = force.select_atoms('all')
    #select_atom_energy = energy.select_atoms('all')
    u2 = mda.Merge(select_atom)
    #loading values from universe
    coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(), select_atom).run().results['timeseries']
    forces = AnalysisFromFunction(lambda ag: ag.positions.copy(), select_atom_force).run().results['timeseries']
    dimensions = AnalysisFromFunction(lambda ag: ag.dimensions.copy(), select_atom).run().results['timeseries']
    #merging it into one the energy is loaded into velocity since it is eeasier to selecting the value by frame/atom and there isn't a way to load custom data into a universe 
    u2.load_new(coordinates, format=MemoryReader, forces=forces, dimensions=dimensions)
    return u2

def main():
    start = dt.datetime.now()
    u2 = load_data() # this should drop all the intermediate values loaded which might benefit memory
    load_data_time = dt.datetime.now()
    print(f"finished loading data: this step = {load_data_time - start}; total ={load_data_time - start}")
    poseidon_object = Poseidon(container=u2, start=1, end=20)
    populate_object_time = dt.datetime.now()
    print(f"finished populate object: this step = {populate_object_time - load_data_time}; total ={populate_object_time - start}")
    analysis_time = dt.datetime.now()
    print(f"finished analysis: this step = {analysis_time - populate_object_time}; total ={analysis_time - start}")
    result = poseidon_object.run_analysis(level_list = ['atomLevel'], verbose=False)
    print("------------------------------SolventData------------------------------")
    print(result[0])
    print("------------------------------SoluteData------------------------------")
    print(result[1])
    try:
        print("------------------------------Contact matrix------------------------------")
        print(result[2])
    except:
        pass

main()