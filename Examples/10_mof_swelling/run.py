import os
from pysimm.apps import mc_md
from pysimm import system

def run(test=False):
    frame = system.read_lammps(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'irmof1_drei.lmps'))
    frame.forcefield = 'dreiding'
    gas1 = system.read_lammps(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ch4.lmps'))
    gas1.forcefield = 'trappe/amber'
    
    mc_props = {'rigid_type': False,
                'max_ins': 1000,
                'Chemical_Potential_Info': -30.09,
                'Temperature_Info': 300,
                'Run_Type': {'steps': 250},
                'CBMC_Info': {'rcut_cbmc': 2.0},
                'Simulation_Length_Info': {'run': 10000,
                                           'coord_freq': 10000,
                                           'prop_freq': 500},
                'VDW_Style': {'cut_val': 9.0},
                'Charge_Style': {'cut_val': 9.0},
                'Property_Info': {'prop1': 'energy_total',
                                  'prop2': 'pressure'}}
    
    md_props = {'temp': 300,
                'pressure': {'start': 15,
                             'iso': 'iso'},
                'timestep': 0.25,
                'cutoff': 9.0,
                'length': 10000,
                'thermo': 2500,
                'dump': 2500,
                'print_to_screen': False}
    
    sim_result = mc_md.mc_md(gas1, frame, mcmd_niter=3, sim_folder='results', mc_props=mc_props, md_props=md_props)
    sim_result.write_lammps('MOFplusME.lmps')
    sim_result.write_xyz('MOFplusME.xyz')

if __name__ == '__main__':
    run(False)