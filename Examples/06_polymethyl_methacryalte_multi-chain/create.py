from pysimm import system, lmps, forcefield, gasteiger
from pysimm.apps.random_walk import random_walk

pmma = system.read_pubchem_smiles('cc(C)(C(=O)OC)')

f = forcefield.Gaff()

pmma.particles[3].linker='head'
pmma.particles[6].linker='tail'

pmma.apply_forcefield(f, charges='gasteiger')

lmps.quick_md(pmma, ensemble='nve')
lmps.quick_min(pmma, min_style='sd')
lmps.quick_min(pmma, min_style='cg')

pmma.write_yaml('pmma_monomer.yaml')

polymer = random_walk(pmma , nmon=5, forcefield=f, density=0.3/4)
polymer = random_walk(pmma , nmon=5, s_=polymer, forcefield=f)
polymer = random_walk(pmma , nmon=5, s_=polymer, forcefield=f)
polymer = random_walk(pmma , nmon=5, s_=polymer, forcefield=f)

polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')
