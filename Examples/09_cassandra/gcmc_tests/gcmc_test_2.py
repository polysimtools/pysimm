from pysimm import system, lmps, cassandra
from collections import OrderedDict
import os

# First, for CASSANDRA GCMC the governing CASSANDRA object should be created
css = cassandra.Cassandra()

# Read the some pre-prepared CASSANDRA .inp parameters file that can be modified then programatically
gcmcPropsRead = css.read_input('my_props.inp')

# Creating fixed polymer system
fixed_sst = system.read_lammps('pim.lmps')
fixed_sst.wrap()

# Setting the one-molecule gas system
gas = system.read_lammps('m-xylene.lmps')

# I Perliminary MC step (pure insertion)
gas_mc = cassandra.McSystem(gas, max_ins=1000, chem_pot=-60.7944)

out_folder = 'gcmc_test_2'
gcmcPropsRead['Run_Name'] = '1.gcmc'
gcmc = cassandra.GCMC(gas_mc, fixed_sst, out_folder=out_folder, **gcmcPropsRead)
css.add_gcmc(gcmc)
css.run()

# old_gcmc.tot_sst.write_lammps('gcmc_test_2/MC-out.lmps')
# gcmc.tot_sst.write_xyz('MC-out.xyz')
# old_gcmc.tot_sst = system.read_lammps('long-MC-out.lmps')

l = 1
while l < 6:
    # >>> 2N: MD step: minimization-dynamics
    sim = lmps.Simulation(gcmc.tot_sst, log=os.path.join(out_folder, str(l) + '.md.log'),
                          print_to_screen=True, cutoff=14.0, dump_unwraped=True)
    sim.add_custom('neighbor 1.0 nsq \nneigh_modify once no every 1 delay 0 check yes')
    fake_md = lmps.MolecularDynamics(ensemble='npt', timestep=0.1, length=4000, thermo=100,
                                     dump=100, dump_name=os.path.join(out_folder, str(l) + '.md.dump'))
    fake_md.write(sim)
    sim.add_custom(gcmc.get_grouped_md(fake_md.input))
    sim.run()

    # Updating the size of the fixed system from the MD simulations for the next MC
    fixed_sst.dim = sim.system.dim

    # >>> 2N + 1: MC step
    new_xyz_file = 'MD' + str(l + 1) + '_out.xyz'
    gcmcPropsRead['Run_Name'] = str(l + 1) + '.gcmc'

    gas_mc = cassandra.McSystem(gas, max_ins=1000, chem_pot=-60.7944)
    # The 'Start_Type' property should be changed in order to make CASSANDRA read previously simulated system
    gcmcPropsRead['Start_Type'] = OrderedDict([('start_type', 'read_config'),
                                               ('species', [1, len(sim.system.molecules) - len(fixed_sst.molecules)]),
                                               ('file_name', os.path.join(out_folder, new_xyz_file))])
    gcmc = cassandra.GCMC(gas_mc, fixed_sst, out_folder=out_folder, **gcmcPropsRead)
    # now we make new .xyz file to read start configuration from
    sim.system.write_xyz(os.path.join(out_folder, new_xyz_file))
    css.add_gcmc(gcmc, is_new=True)
    css.run()
    l += 1
