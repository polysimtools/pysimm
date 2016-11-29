from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import copolymer

# use a smiles string to query the pubchem search database and read the mol file returned from the http request
pe = system.read_pubchem_smiles('cc')

# we'll instantiate a GAFF2 forcefield object for use later
f = forcefield.Gaff2()

# particles 1 and 2 in the monomer are going to be the head and tail linkers
pe.particles[1].linker='head'
pe.particles[2].linker='tail'

# the resulting system has sufficient information to type with the forcefield object we made earlier
# we will also determine partial charges using the gasteiger algorithm
pe.apply_forcefield(f, charges='gasteiger')

# do a quick minimization of the monomer
lmps.quick_min(pe, min_style='fire')

# write a yaml file for the pe monomer
pe.write_yaml('pe_monomer.yaml')

# use a smiles string to query the pubchem search database and read the mol file returned from the http request
ps = system.read_pubchem_smiles('cc(C1=CC=CC=C1)')

# particles 8 and 7 in the monomer are going to be the head and tail linkers
ps.particles[8].linker='head'
ps.particles[7].linker='tail'

# like in example 3, we need to identify the bonds in the ring as aromatic
for b in ps.bonds:
  if (not b.a.linker and not b.b.linker) and b.a.elem=='C' and b.b.elem=='C':
    b.order='A'

# the resulting system has sufficient information to type with the forcefield object we made earlier
# we will also determine partial charges using the gasteiger algorithm
ps.apply_forcefield(f, charges='gasteiger')

# do a quick minimization of the monomer
lmps.quick_min(ps, min_style='fire')

# write a yaml file for the ps monomer
ps.write_yaml('ps_monomer.yaml')

# run the copolymer random walk method with 10 total repeat units, using an alternating pattern
# here, settings is passed to the LAMMPS simulations during the polymerization, and we will use 2 processors
polymer = copolymer([pe, ps], 10, pattern=[1, 1], forcefield=f, settings={'np': 2})

# write a few different file formats
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')