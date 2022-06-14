from pysimm import system, lmps, forcefield


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
    
    molecule_list = [ethanol, acetone]
    
    if test:
        n_molecules = [20, 20]
    else:
        n_molecules = [200, 200]
    
    s = system.replicate(molecule_list, n_molecules, density=0.3)
    
    min_settings = {
        'name': 'cg_min',
        'min_style': 'cg',
        'maxiter': int(5e+5),
        'maxeval': int(5e+6),
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
        'length': 2500
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
        'length': 5000,
    }
    
    npt_settings_add = {
        'name': 'npt_md',
        'print_to_screen': True,
        'ensemble': 'npt',
        'temperature': 300,
        'new_v': True,
        'pressure': {
            'start': 1,
            'stop': 1
        },
        'length': 5000,
    }
    
    if test:
        nvt_settings['length'] = 2000
        npt_settings['length'] = 2000
        
    sim = lmps.Simulation(s)
    sim.add_min(**min_settings)
    
    sim.add(lmps.OutputSettings(thermo={'freq': 500, 
                                        'style': 'custom', 
                                        'args': ['step', 'temp', 'etotal', 'press', 'density']}))

    sim.add_md(**nvt_settings)
    sim.add_md(**npt_settings)
    sim.add_md(**npt_settings_add)
    sim.run()

    s.write_yaml('mixture.yaml')
    s.write_lammps('mixture.lmps')


if __name__ == '__main__':
    run()
