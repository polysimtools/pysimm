from pysimm import system, lmps, forcefield, gasteiger
from pysimm.apps.random_walk import random_walk

# use a smiles string to query the pubchem search database and read the mol file returned from the http request
pmma = system.read_pubchem_smiles('cc(C)(C(=O)OC)')

# we'll instantiate a GAFF forcefield object for use later
f = forcefield.Gaff()

# particles 3 and 6 in the monomer are going to be the head and tail linkers
pmma.particles[3].linker='head'
pmma.particles[6].linker='tail'

# the resulting system has sufficient information to type with the forcefield object we made earlier
# we will also determine partial charges using the gasteiger algorithm
pmma.apply_forcefield(f, charges='gasteiger')

# do a quick minimization of the monomer
lmps.quick_min(pmma, min_style='fire')

# write a yaml file for the pmma monomer
pmma.write_yaml('pmma_monomer.yaml')

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
