from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk

# use a smiles string to query the pubchem search database and read the mol file returned from the http request
pmma = system.read_pubchem_smiles('cc(C)(C(=O)OC)')

# particles 3 and 6 in the monomer are going to be the head and tail linkers
pmma.particles[3].linker='head'
pmma.particles[6].linker='tail'

# we'll instantiate a Dreiding forcefield object for use later
f = forcefield.Dreiding()

# the resulting system has sufficient information to type with the forcefield object we made earlier
# we will also determine partial charges using the gasteiger algorithm
pmma.apply_forcefield(f, charges='gasteiger')

# do a quick minimization of the monomer
lmps.quick_min(pmma, min_style='fire')

pmma.pair_style = 'lj'

# we're going to make 4 chains, each of 5 repeat units
# the first system we make will be used as the initial system for the subsequent random walk calls
polymer = random_walk(pmma , nmon=5, forcefield=f, density=0.3/4)
polymer = random_walk(pmma , nmon=5, s_=polymer, forcefield=f)
polymer = random_walk(pmma , nmon=5, s_=polymer, forcefield=f)
polymer = random_walk(pmma , nmon=5, s_=polymer, forcefield=f)

# write a few different file formats
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')
