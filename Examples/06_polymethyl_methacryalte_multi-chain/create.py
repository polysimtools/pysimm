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
pmma.write_lammps('pmma_monomer.lmps')

pmma.viz()

polymer = random_walk(pmma , nmon=10, forcefield=f, density=0.3/4, settings={'np': 2})
polymer = random_walk(pmma , nmon=10, s_=polymer, forcefield=f, settings={'np': 4})
polymer = random_walk(pmma , nmon=10, s_=polymer, forcefield=f, settings={'np': 6})
polymer = random_walk(pmma , nmon=10, s_=polymer, forcefield=f, settings={'np': 8})

polymer.viz()
