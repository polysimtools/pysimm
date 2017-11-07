from pysimm import system, lmps, forcefield
# from pysimm import amber

def run(test=False):
    try:
        ethanol = system.read_pubchem_smiles('CCO')
    except:
        import os
        ethanol = system.read_mol(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'CCO.mol'))
    try:
        acetone = system.read_pubchem_smiles('CC(=O)C')
    except:
        import os
        acetone = system.read_mol(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'CC(=O)C.mol'))
    
    f = forcefield.Gaff2()
    
    ethanol.apply_forcefield(f, charges='gasteiger')
    acetone.apply_forcefield(f, charges='gasteiger')
    
    # amber.calc_charges(ethanol)
    # amber.calc_charges(acetone)
    
    lmps.quick_min(ethanol, min_style='fire')
    lmps.quick_min(acetone, min_style='fire')
    
    molecule_list=[ethanol,acetone]
    
    if test:
        n_molecules=[30,20]
    else:
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
        'temperature': {
            'start': 100,
            'stop': 300
        },
        'new_v': True,
        'length': 10000
    }
    
    npt_settings = {
        'name': 'npt_md',
        'print_to_screen': True,
        'ensemble': 'npt',
        'temperature': 300,
        'new_v': True,
        'pressure': {
            'start': 1000,
            'stop': 1
        },
        'length': 100000,
        'thermo_style': 'custom step temp press density'
    }
    
    if test:
        nvt_settings['length'] = 2000
        npt_settings['length'] = 2000
        
    
    lmps.quick_min(s, **min_settings)
    lmps.quick_md(s, **nvt_settings)
    lmps.quick_md(s, **npt_settings)
    
    s.write_xyz('mixture.xyz')
    s.write_yaml('mixture.yaml')
    s.write_lammps('mixture.lmps')
    s.write_chemdoodle_json('mixture.json')
    
if __name__ == '__main__':
    run()