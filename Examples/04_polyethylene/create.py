from pysimm import system, lmps, forcefield, gasteiger
from pysimm.apps.random_walk import random_walk

s = system.read_pubchem_smiles('cc')

f = forcefield.Gaff()

s.particles[1].linker='head'
s.particles[2].linker='tail'

s.apply_forcefield(f, charges='gasteiger')

lmps.quick_min(s, min_style='fire')

s.write_xyz('pe_mnomer.xyz')
s.write_yaml('pe_monomer.yaml')
s.write_lammps('pe_monomer.lmps')
s.write_chemdoodle_json('pe_monomer.json')

polymer = random_walk(s, 10, forcefield=f)

polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')