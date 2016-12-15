from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk
from pysimm.models.dreiding.pe import monomer

# we'll create a pe monomer from the pysimm.models database
pe = monomer()

# we'll instantiate a Dreiding forcefield object for use later
f = forcefield.Dreiding()

# run the random_walk polymerization method making a chain of 10 repeat units
# the forcefield object is supplied to get new forcefield types not in the monomer system
polymer = random_walk(pe, 10, forcefield=f)

# write a few different file formats
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')
