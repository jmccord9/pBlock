import numpy as np
import os
from ase.units import kJ, mol
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.io.trajectory import TrajectoryWriter
from ase.visualize import view
from ase.io import read, write
from ase.build import add_adsorbate
from ase.constraints import Hookean
from ase.constraints import FixAtoms

def read_poscar(path):
    slab = read(path)
    # automatically fix all atoms when reading in a POSCAR
    # fix = FixAtoms(indices = [atom.index for atom in sim if atom.symbol == 'Ti' or atom.symbol == 'O'])
    # slab.set_constraint(fix)
    return slab

def write_obj(write_path,obj):
    # write .traj or POSCAR files
    write(write_path,obj)
    return None

def place_dopant(atoms,dopantAdd,dopantThere):
    symbols = "Ti2O4Ti2O4Ti2O4Ti2O4Ti2O4Ti2O4Ti2O4Ti2O"+dopantAdd+"O2"
    pos, con, cell = atoms.positions, atoms.constraints, atoms.cell
    new_atoms = Atoms(symbols=symbols,pbc=True,cell=cell,positions=pos,constraint=con)
    return new_atoms

def fix_atoms(atoms,constraint):
    # fix = FixAtoms(indices = [atom.index for atom in sim if atom.symbol == 'Ti' or atom.symbol == 'O'])
    # slab.set_constraint(fix)
    fix = FixAtoms(mask=constraint)
    atoms.set_constraint(fix)
    return atoms

def place_adsorbate(slab,species):
    add_adsorbate(slab, species, 1.34823, (-0.2305603365560289, 2.7825235214147974))
    add_adsorbate(slab, species, 2.48744, (-0.3525537826397047, 2.6017940132659629))
    return None

def view_atoms(atoms_obj):
    view(atoms_obj)
    return None

def get_pe(atoms_obj):
    return atoms_obj.get_potential_energy()

def read_traj(traj_file):
    return read(traj_file, index=-1)

# fixDope = FixAtoms(indices = [atom.index for atom in sim])
# hook = Hookean(a1 = 48, a2 = 49, rt = 1.09, k = 14)
# sim.set_constraint(hook)
# sim.set_constraint(fixDope)

if __name__ == '__main__':

    base = '/Users/jamesmccord/Dropbox (GaTech)/pBlock'
    sem_dir = 'fall2021'
    sub_dir = 'Br'

    read_path = os.path.join(base,sem_dir,'POSCAR')
    write_path = os.path.join(base,sem_dir,'slab_base.traj')

    slab = read_poscar(read_path)
    slab = fix_atoms(slab,slab.positions[:, 2] < 11.0)
    slab = place_dopant(slab,'Br','B')
    write_obj(write_path, slab)
    view_atoms(slab)

