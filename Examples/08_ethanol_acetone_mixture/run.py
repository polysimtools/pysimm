from pysimm import system, lmps, forcefield
# from pysimm import amber

ethanol = system.read_pubchem_smiles('CCO')
acetone = system.read_pubchem_smiles('CC(=O)C')

f = forcefield.Gaff2()

ethanol.apply_forcefield(f, charges='gasteiger')
acetone.apply_forcefield(f, charges='gasteiger')

# amber.calc_charges(ethanol)
# amber.calc_charges(acetone)

lmps.quick_min(ethanol, min_style='fire')
lmps.quick_min(acetone, min_style='fire')

molecule_list=[ethanol,acetone]
n_molecules=[300,200]

s=system.replicate(molecule_list, n_molecules , density=0.3)

lmps.quick_min(s, name='fire_min', min_style='fire', print_to_screen=True, thermo=500)
lmps.quick_md(s, name='nvt_md', ensemble='nvt', length=10000, print_to_screen=True, thermo=500)
lmps.quick_md(s, name='npt_md', ensemble='npt', length=100000, print_to_screen=True, thermo=500)

s.write_xyz('mixture.xyz')
s.write_yaml('mixture.yaml')
s.write_lammps('mixture.lmps')
s.write_chemdoodle_json('mixture.json')