from pysimm import system, lmps, cassandra
from collections import OrderedDict
import os
import random

MAX_INS = [2000, 2000]
CHEM_POT = [-33.15, -60.37]
TMPR = 300
L_MAX = 6

# First, for CASSANDRA GCMC the governing CASSANDRA object should be created
css = cassandra.Cassandra()

# Read the some pre-prepared CASSANDRA .inp parameters file that can be modified then programatically
gcmcPropsRead = css.read_input('my_props.inp')

# Creating fixed polymer system
fixed_sst = system.read_lammps('pim.lmps')
fixed_sst.wrap()

# Setting the one-molecule gas system
co2 = system.read_lammps('co2.lmps')
xylene = system.read_lammps('m-xylene.lmps')

# I Preliminary MC step (pure insertion)
gas_mc = cassandra.McSystem([co2, xylene], chem_pot=CHEM_POT, max_ins=MAX_INS)
out_folder = 'gcmc_test_2'
gcmcPropsRead['Run_Name'] = '1.gcmc'
gcmc = cassandra.GCMC(gas_mc, fixed_sst, out_folder=out_folder, props_file='1.gcmc_props.inp', **gcmcPropsRead)
css.add_gcmc(gcmc)
css.run()

l = 1
time_steps = [0.15] * L_MAX
while l < L_MAX:
    # >>> 2N: MD (LAMMPS) step:
    sim_sst = gcmc.tot_sst
    sim_sst.write_lammps(os.path.join(out_folder, str(l) + '.before_md.lmps'))
    sim = lmps.Simulation(sim_sst, log=os.path.join(out_folder, str(l) + '.md.log'),
                          print_to_screen=True, cutoff=14.0)

    # custom definitions for the neighbour list updates
    sim.add_custom('neighbor 1.0 nsq \nneigh_modify once no every 1 delay 0 check yes')

    # custom definition of md with rigid gas molecules
    tmp_md = lmps.MolecularDynamics(ensemble='npt', timestep=time_steps[l - 1], length=10000, thermo=100,
                                    temp=TMPR, dump=100, dump_name=os.path.join(out_folder, str(l) + '.md.dump'))
    tmp_md.write(sim)
    custom_inp = gcmc.get_grouped_md(tmp_md.input)

    # custom simulation input for velocities restriction
    line = 'velocity {:} create {:} {:}\n'
    vel_def = line.format(gcmc.fxd_sst.name, TMPR, random.randint(int(1e+6), int(1e+7 - 1)))
    vel_def += line.format(gcmc.mc_sst.name, TMPR, random.randint(int(1e+6), int(1e+7 - 1)))
    vel_def += 'run 0\n'
    vel_def += line.replace('create', 'scale').format(gcmc.fxd_sst.name, TMPR, '')
    vel_def += line.replace('create', 'scale').format(gcmc.mc_sst.name, TMPR, '')

    # adding all customs to LAMMPS input and run simulations
    sim.add_custom(custom_inp.replace('\nrun', '\n' + vel_def + 'run'))
    sim.run()

    # Updating the size of the fixed system from the MD simulations for the next MC
    fixed_sst.dim = sim.system.dim

    # >>> 2N + 1: MC (CASSANDRA) step
    new_xyz_file = 'MD' + str(l) + '_out.xyz'
    gcmcPropsRead['Run_Name'] = str(l + 1) + '.gcmc'

    # The 'Start_Type' property should be changed in order to make CASSANDRA read previously simulated system
    gcmcPropsRead['Start_Type'] = OrderedDict([('start_type', 'read_config'), ('species', [1] + gas_mc.made_ins),
                                               ('file_name', os.path.join(out_folder, new_xyz_file))])

    gas_mc = cassandra.McSystem([co2, xylene], chem_pot=CHEM_POT, max_ins=MAX_INS)
    gcmc = cassandra.GCMC(gas_mc, fixed_sst, out_folder=out_folder, props_file=str(l + 1) + '.gcmc_props.inp',
                          **gcmcPropsRead)
    # now we make new .xyz file to read start configuration from
    sim.system.write_xyz(os.path.join(out_folder, new_xyz_file))
    css.add_gcmc(gcmc, is_new=True)
    css.run()
    l += 1
