from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk
from pysimm.models.monomers.gaff2.pe import monomer

# we'll create a pe monomer from the pysimm.models database
pe = monomer()

# we'll instantiate a GAFF2 forcefield object for use later
f = forcefield.Gaff2()

# the monomers do not have any charges, so we will derive partial charges using the gasteiger algorithm
pe.apply_charges(f, charges='gasteiger')

# run the random_walk polymerization method making a chain of 10 repeat units
# the forcefield object is supplied to get new forcefield types not in the monomer system
polymer = random_walk(pe, 10, forcefield=f)

# write a few different file formats
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')

# if you want to restart a polymerization, the yaml file format retains linker information
# random_walk looks for the last head and tail linkers, so just run a copolymerization with the original polymer chain and new monomers
# we give the copolymer function a list of reference "monomers", but use the first polymer chain as the first "monomer" and only insert one
# then we use the pattern argument to define how many of each "monomers" to add. Let's add 5 more monomers to our chain

# first import the copolymer function

from pysimm.apps.random_walk import copolymer

# now read in the yaml file we saved after making our first polymer

original_polymer = system.read_yaml('polymer.yaml')

# we can use our original polyethylene monomer because it doesn't get modified during polymerization
# the total number of monomers we're adding is 6, 1 for the original polymer chain, and 5 for our new monomers

longer_polymer = copolymer([original_polymer, pe], 6, pattern=[1, 5], forcefield=f)

longer_polymer.write_xyz('longer_polymer.xyz')
longer_polymer.write_yaml('longer_polymer.yaml')
longer_polymer.write_lammps('longer_polymer.lmps')
longer_polymer.write_chemdoodle_json('longer_polymer.json')