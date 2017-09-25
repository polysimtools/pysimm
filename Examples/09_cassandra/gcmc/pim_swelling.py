from pysimm.apps import mc_md
from pysimm import system, cassandra

polymer = system.read_lammps('pim.lmps')

cs = cassandra.Cassandra(polymer)
my_mc_props = cs.read_input('my_props.inp')
my_mc_props.update({'rigid_type': [True, False],
                    'max_ins': [2000, 1500],
                    'Chemical_Potential_Info': [-35.01, -26.15]})

my_md_props = {'ensemble': 'npt',
               'timestep': 0.5,
               'length': 1e+4,
               'thermo': 2.5e+3,
               'temp': 300,
               'dump': 2.5e+3,
               'print_to_screen': False,
               'cutoff': 14.0}

gas1 = system.read_lammps('co2.lmps')
gas2 = system.read_lammps('ch4.lmps')

mc_md.mc_md([gas1, gas2], polymer,
            mcmd_niter=21,
            sim_folder='results',
            mc_props=my_mc_props,
            md_props=my_md_props)
