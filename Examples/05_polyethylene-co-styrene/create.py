from pysimm import system, lmps, forcefield, gasteiger
from pysimm.apps.random_walk import copolymer

pe = system.read_pubchem_smiles('cc')

f = forcefield.Gaff()

pe.particles[1].linker='head'
pe.particles[2].linker='tail'

pe.apply_forcefield(f, charges='gasteiger')

lmps.quick_md(pe, ensemble='nve')
lmps.quick_min(pe, min_style='sd')
lmps.quick_min(pe, min_style='cg')

pe.write_yaml('pe_monomer.yaml')
pe.write_lammps('pe_monomer.lmps')

pe.viz()

ps = system.read_pubchem_smiles('cc(C1=CC=CC=C1)')

ps.particles[8].linker='head'
ps.particles[7].linker='tail'

for b in ps.bonds:
  if (not b.a.linker and not b.b.linker) and b.a.elem=='C' and b.b.elem=='C':
    b.order='A'

ps.apply_forcefield(f, charges='gasteiger')

lmps.quick_md(ps, ensemble='nve')
lmps.quick_min(ps, min_style='sd')
lmps.quick_min(ps, min_style='cg')

ps.write_yaml('ps_monomer.yaml')
ps.write_lammps('ps_monomer.lmps')

ps.viz()

polymer = copolymer([pe, ps], 10, pattern=[1, 1], forcefield=f, settings={'np': 2})

polymer.viz()
