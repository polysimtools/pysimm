from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk

# use a smiles string to query the pubchem search database and read the mol file returned from the http request
s = system.read_pubchem_smiles('cc')

# we'll instantiate a Dreiding forcefield object for use later
f = forcefield.Dreiding()

# particles 1 and 2 in the monomer are going to be the head and tail linkers
s.particles[1].linker='head'
s.particles[2].linker='tail'

# the resulting system has sufficient information to type with the forcefield object we made earlier
# we will also determine partial charges using the gasteiger algorithm
s.apply_forcefield(f, charges='gasteiger')

# do a quick minimization of the monomer
lmps.quick_min(s, min_style='fire')

# write a few different file formats
s.write_xyz('pe_monomer.xyz')
s.write_yaml('pe_monomer.yaml')
s.write_lammps('pe_monomer.lmps')
s.write_chemdoodle_json('pe_monomer.json')

# run the random_walk polymerization method making a chain of 10 repeat units
# the forcefield object is supplied to get new forcefield types not in the monomer system
polymer = random_walk(s, 10, forcefield=f)

# write a few different file formats
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')
