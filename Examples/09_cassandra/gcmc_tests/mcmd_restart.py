from pysimm import system, lmps, cassandra
from collections import OrderedDict
import os
import re
import random

MAX_INS = 2000
CHEM_POT = -33.15
TMPR = 300
L_MAX = 11

# Restarting parameters
restart_ind = 11

# First, for CASSANDRA GCMC the governing CASSANDRA object should be created
css = cassandra.Cassandra()
# Read the some pre-prepared CASSANDRA .inp parameters file that can be modified then programatically
gcmcPropsRead = css.read_input('my_props.inp')
# Creating fixed polymer system
fixed_sst = system.read_lammps('pim.lmps')
# Reading "reastart dump" that is simply the .lmps file from one of the previous MC-MD iterations
fixed_ref_sst = system.read_lammps(os.path.join('results', '10.md.lammps'))
for p in fixed_sst.particles:
    p.x = fixed_ref_sst.particles[p.tag].x
    p.y = fixed_ref_sst.particles[p.tag].y
    p.z = fixed_ref_sst.particles[p.tag].z
fixed_sst.dim = fixed_ref_sst.dim

# Setting the gas system
co2 = system.read_lammps('co2.lmps')

# I Reading and arraging data from previous MC step
gas_mc = cassandra.McSystem(co2, chem_pot=CHEM_POT, max_ins=MAX_INS, is_rigid=True)
out_folder = 'results'
gcmcPropsRead['Run_Name'] = str(restart_ind) + '.gcmc'
gcmc = cassandra.GCMC(gas_mc, fixed_sst, out_folder=out_folder, props_file=str(restart_ind) + '.gcmc_props.inp',
                      **gcmcPropsRead)
gcmc.upd_simulation()

l = restart_ind  # loop variable
time_step = 0.5
lngth = 2e+6
nonrig_group_name = 'nonrigid_b'
rig_group_name = 'rigid_b'
while l < L_MAX + restart_ind:
    # >>> 2N: MD (LAMMPS) step:
    sim = lmps.Simulation(gcmc.tot_sst, log=os.path.join(out_folder, str(l) + '.md.log'),
                          print_to_screen=True, cutoff=14.0)

    # custom definitions for the neighbour list updates
    sim.add_custom('neighbor 1.0 nsq \nneigh_modify once no every 1 delay 0 check yes')

    # adding group definitions to separate rigid and non-rigid bodies
    grp_tmpl = 'group {:} id {:}'
    sim.add_custom(grp_tmpl.format(nonrig_group_name, gcmc.group_by_id('nonrigid')[0]))
    sim.add_custom(grp_tmpl.format(rig_group_name, gcmc.group_by_id('rigid')[0]))

    # create the description of the molecular dynamics simulation
    tmp_md = lmps.MolecularDynamics(ensemble='npt', timestep=time_step, length=int(lngth), thermo=2500, pressure=20,
                                    temp=TMPR, dump=5000, dump_name=os.path.join(out_folder, str(l) + '.md.dump'),
                                    new_v=True, scale_v=True)

    # obtain the simulation (LAMMPS) input in order to customly modify it later
    tmp_md.write(sim)

    # replace single default fix with two separate fixes for rigid and nonrigid bodies separately
    old_line = re.search('(?<=(\nfix)).*', tmp_md.input).group(0)
    corr_fix = re.sub('all', nonrig_group_name, old_line) + ' dilate all\n'
    corr_fix += 'fix' + re.sub('iso\s+\d+[.\d+]*\s+\d+[.\d+]*\s+\d+[.\d+]*', '', old_line). \
        replace('1', '2', 1). \
        replace('all', rig_group_name). \
        replace('npt', 'rigid/nvt/small molecule') + '\n'

    # adding the spring fix to the geometrical center of the system to avoid system creep
    corr_fix += 'fix {:} {:} spring tether {:} {:} {:} {:} {:}\n'.format(3, nonrig_group_name, 30.0, 0.0, 0.0, 0.0, 0.0)

    # saving all fixes to the input
    tmp_md.input = tmp_md.input.replace(old_line, corr_fix)

    # adding "run 0" line for correct temperature scaling of the system with rigid molecules
    tmp_md.input = tmp_md.input.replace('velocity all scale', 'run 0\nvelocity all scale')

    # The input for correct simulations is set, starting LAMMPS:
    sim.add_custom(tmp_md.input)
    sim.run()

    # Updating the size of the fixed system from the MD simulations for the next MC
    fixed_sst.dim = sim.system.dim

    # >>> 2N + 1: MC (CASSANDRA) step
    new_xyz_file = 'MD' + str(l) + '_out.xyz'
    gcmcPropsRead['Run_Name'] = str(l + 1) + '.gcmc'

    # The 'Start_Type' property should be changed in order to make CASSANDRA read previously simulated system
    gcmcPropsRead['Start_Type'] = OrderedDict([('start_type', 'read_config'),
                                               ('species', [1] + [0] * len(gas_mc.chem_pot)),
                                               ('file_name', os.path.join(out_folder, new_xyz_file))])

    gas_mc = cassandra.McSystem(co2, chem_pot=CHEM_POT, max_ins=MAX_INS, is_rigid=True)
    gcmc = cassandra.GCMC(gas_mc, fixed_sst, out_folder=out_folder,
                          props_file=str(l + 1) + '.gcmc_props.inp', **gcmcPropsRead)
    # now we make new .xyz file to read start configuration from
    sim.system.write_xyz(os.path.join(out_folder, new_xyz_file))
    css.add_gcmc(gcmc, is_new=True)
    css.run()
    sim_sst = gcmc.tot_sst
    l += 1
