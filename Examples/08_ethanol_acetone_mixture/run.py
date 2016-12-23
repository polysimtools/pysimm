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

min_settings = {
    'name': 'fire_min',
    'min_style': 'fire',
    'print_to_screen': True
}

nvt_settings = {
    'name': 'nvt_md',
    'print_to_screen': True,
    'ensemble': 'nvt',
    't_start': 100,
    't_stop': 300,
    'new_v': True,
    'length': 10000
}

npt_settings = {
    'name': 'npt_md',
    'print_to_screen': True,
    'ensemble': 'npt',
    'temp': 300,
    'new_v': True,
    'p_start': 1000,
    'p_stop': 1,
    'length': 100000,
    'thermo_style': 'custom step temp press density'
}

lmps.quick_min(s, **min_settings)
lmps.quick_md(s, **nvt_settings)
lmps.quick_md(s, **npt_settings)

s.write_xyz('mixture.xyz')
s.write_yaml('mixture.yaml')
s.write_lammps('mixture.lmps')
s.write_chemdoodle_json('mixture.json')