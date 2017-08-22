import re

from pysimm import system, lmps, cassandra
from collections import OrderedDict
import os

MAX_INS = [2000, 2000]
CHEM_POT = [-33.15, -60.37]
TMPR = 300
L_MAX = 11

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
gas_mc = cassandra.McSystem([co2, xylene], chem_pot=CHEM_POT, max_ins=MAX_INS, is_rigid=[True, False])
out_folder = 'gcmc_test_2'
gcmcPropsRead['Run_Name'] = '1.gcmc'
gcmc = cassandra.GCMC(gas_mc, fixed_sst, out_folder=out_folder, props_file='1.gcmc_props.inp', **gcmcPropsRead)
css.add_gcmc(gcmc)
css.run()

l = 1
time_step = 0.5
lngth = 2e+4
nonrig_group_name = 'nonrigid_b'
rig_group_name = 'rigid_b'

while l < L_MAX:
    # >>> 2N: MD (LAMMPS) step:
    sim_sst = gcmc.tot_sst
    sim_sst.write_lammps(os.path.join(out_folder, str(l) + '.before_md.lmps'))
    sim = lmps.Simulation(sim_sst, log=os.path.join(out_folder, str(l) + '.md.log'),
                          print_to_screen=True, cutoff=14.0)

    # custom definitions for the neighbour list updates
    sim.add_custom('neighbor 1.0 nsq \nneigh_modify once no every 1 delay 0 check yes')

    # adding group definitions to separate rigid and non-rigid bodies
    grp_tmpl = 'group {:} id {:}'
    sim.add_custom(grp_tmpl.format(nonrig_group_name, gcmc.nonrigid_idxs))
    sim.add_custom(grp_tmpl.format(rig_group_name, gcmc.rigid_idxs))

    # create the description of the molecular dynamics simulation
    tmp_md = lmps.MolecularDynamics(ensemble='npt', timestep=time_step, length=int(lngth), thermo=1000,
                                    temp=TMPR, dump=1000, dump_name=os.path.join(out_folder, str(l) + '.md.dump'),
                                    new_v=True, scale_v=True)

    # obtain the simulation (LAMMPS) input in order to customly modify it later
    tmp_md.write(sim)

    # replace single default fix with two separate fixes for rigid and nonrigid bodies separately
    old_line = re.search('(?<=(\nfix)).*', tmp_md.input).group(0)
    corr_fix = re.sub('all', nonrig_group_name, old_line) + ' dilate all\n'

    corr_fix += 'fix' + re.sub('iso\s+\d+.\d+\s+\d+.\d+\s+\d+.\d+', '', old_line).\
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
    sim.run(np=8)

    # Updating the size of the fixed system from the MD simulations for the next MC
    fixed_sst.dim = sim.system.dim

    # >>> 2N + 1: MC (CASSANDRA) step
    new_xyz_file = 'MD' + str(l) + '_out.xyz'
    gcmcPropsRead['Run_Name'] = str(l + 1) + '.gcmc'

    # The 'Start_Type' property should be changed in order to make CASSANDRA read previously simulated system
    gcmcPropsRead['Start_Type'] = OrderedDict([('start_type', 'read_config'),
                                               ('species', [1] + [0] * len(gas_mc.chem_pot)),
                                               ('file_name', os.path.join(out_folder, new_xyz_file))])

    gas_mc = cassandra.McSystem([co2, xylene], chem_pot=CHEM_POT, max_ins=MAX_INS)
    gcmc = cassandra.GCMC(gas_mc, fixed_sst, out_folder=out_folder,
                          props_file=str(l + 1) + '.gcmc_props.inp', **gcmcPropsRead)
    # now we make new .xyz file to read start configuration from
    sim.system.write_xyz(os.path.join(out_folder, new_xyz_file))
    css.add_gcmc(gcmc, is_new=True)
    css.run()
    l += 1
