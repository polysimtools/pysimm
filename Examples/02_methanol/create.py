from pysimm import system, lmps, forcefield, gasteiger

s = system.read_pubchem_smiles('CO')

s.apply_forcefield(forcefield.Gaff())

s.write_yaml('methanol.yaml')
s.write_lammps('methanol.lmps')

lmps.quick_min(s)
lmps.quick_min(s, min_style='cg')

s.viz()
