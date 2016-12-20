from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk
from pysimm.models.monomers.gaff2.pmma import monomer

# we'll create a pmma monomer from the pysimm.models database
pmma = monomer()

# we'll instantiate a GAFF2 forcefield object for use later
f = forcefield.Gaff2()

# we're going to make 4 chains, each of 5 repeat units
# the first system we make will be used as the initial system and then replicated to form 4 chains
# in this case the system.replicate function take a system input(polymer),
# replicates a defined number of times(4), and inserts the new replication randomly(rand=True) at the specified density(0.022) 

print('Building polymer chain 1...')
polymer = random_walk(pmma , nmon=5, forcefield=f)
print('Replicating polymer chain...')
uniform_polymer = system.replicate(polymer, 4, density = 0.022, rand=True)

# next we're going to make 4 chains, with lengths of 2, 4, 6, and 8 monomer units
# the first system we make will be used as the initial system for the subsequent random walk calls

print('Building polymer chain 1...')
nonuniform_polymer = random_walk(pmma , nmon=2, forcefield=f, density=0.3/4)
print('Building polymer chain 2...')
nonuniform_polymer = random_walk(pmma , nmon=4, s_=nonuniform_polymer, forcefield=f)
print('Building polymer chain 3...')
nonuniform_polymer = random_walk(pmma , nmon=6, s_=nonuniform_polymer, forcefield=f)
print('Building polymer chain 4...')
nonuniform_polymer = random_walk(pmma , nmon=8, s_=nonuniform_polymer, forcefield=f)

# now that we have our two polymer systems, let's calculate their molecular weight dispersity
uniform_polymer.set_mm_dist()
nonuniform_polymer.set_mm_dist()

print('')
print('Uniform polymer')
print('---------------')
print('Number average molecular weight: {}'.format(uniform_polymer.mn))
print('Weight average molecular weight: {}'.format(uniform_polymer.mw))
print('Dispersity: {}'.format(uniform_polymer.dispersity))
print('')
print('Nonuniform polymer')
print('------------------')
print('Number average molecular weight: {}'.format(nonuniform_polymer.mn))
print('Weight average molecular weight: {}'.format(nonuniform_polymer.mw))
print('Dispersity: {}'.format(nonuniform_polymer.dispersity))

# write a few different file formats
uniform_polymer.write_yaml('uniform_polymer.yaml')
uniform_polymer.write_xyz('uniform_polymer.xyz')

nonuniform_polymer.write_yaml('nonuniform_polymer.yaml')
nonuniform_polymer.write_xyz('nonuniform_polymer.xyz')
