from pysimm import system, lmps, forcefield, gasteiger

s = system.read_pubchem_smiles('CO')

s.apply_forcefield(forcefield.Gaff())

s.write_xyz('methanol.xyz')
s.write_yaml('methanol.yaml')
s.write_lammps('methanol.lmps')
s.write_chemdoodle_json('methanol.json')

# we'll perform energy minimization using the fire algorithm in LAMMPS
lmps.quick_min(s, min_style='fire')
