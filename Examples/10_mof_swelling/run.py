from pysimm.apps import mc_md
from pysimm import system


def run(test=False):
    frame = system.read_lammps('irmof-14.lmps')
    frame.forcefield = 'dreiding'
    gas1 = system.read_lammps('ch4.lmps')
    gas1.forcefield = 'trappe/amber'
    
    mc_props = {'rigid_type': False,
                'max_ins': 1000,
                'Chemical_Potential_Info': -22.5037,
                'Temperature_Info': 300,
                'Rcutoff_Low': 0.0,
                'Run_Type': {'steps': 100},
                'CBMC_Info': {'rcut_cbmc': 2.0},
                'Simulation_Length_Info': {'run': 300000,
                                           'coord_freq': 100000,
                                           'prop_freq': 1000},
                'VDW_Style': {'cut_val': 14.0},
                'Charge_Style': {'cut_val': 14.0},
                'Property_Info': {'prop1': 'energy_total',
                                  'prop2': 'pressure',
                                  'prop3': 'nmols'}}
    
    md_props = {'temp': 300,
                'pressure': {'start': 300,
                             'iso': 'iso'},
                'timestep': 1,
                'cutoff': 14.0,
                'length': 100000,
                'thermo': 2500,
                'dump': 25000,
                'np': 4,
                'print_to_screen': True}
    
    sim_result = mc_md.mc_md(gas1, frame, mcmd_niter=5, sim_folder='results', mc_props=mc_props, md_props=md_props)
    sim_result.write_lammps('MOFplusME.lmps')
    sim_result.write_xyz('MOFplusME.xyz')

if __name__ == '__main__':
    run(False)

