from pysimm import system, lmps, forcefield, gasteiger

s = system.read_pubchem_smiles('c1=cc=cc=c1')

for b in s.bonds:
  if b.a.elem=='C' and b.b.elem=='C':
    b.order='A'

s.apply_forcefield(forcefield.Gaff(), charges='gasteiger')

s.write_xyz('benzene.xyz')
s.write_yaml('benzene.yaml')
s.write_lammps('benzene.lmps')
s.write_chemdoodle_json('benzene.json')

# we'll perform energy minimization using the fire algorithm in LAMMPS
lmps.quick_min(s, min_style='fire')
