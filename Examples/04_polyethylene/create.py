from pysimm import system, lmps, forcefield, gasteiger
from pysimm.apps.random_walk import random_walk

s = system.read_pubchem_smiles('cc')

f = forcefield.Gaff()

s.particles[1].linker='head'
s.particles[2].linker='tail'

s.apply_forcefield(f, charges='gasteiger')

lmps.quick_md(s, ensemble='nve')
lmps.quick_min(s, min_style='sd')
lmps.quick_min(s, min_style='cg')

s.write_yaml('pe_monomer.yaml')
s.write_lammps('pe_monomer.lmps')

s.viz()

polymer = random_walk(s, 10, forcefield=f)

polymer.viz()
