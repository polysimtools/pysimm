import os
import re
from collections import OrderedDict
from pysimm import system, lmps, cassandra
import glob

MAX_INS = 2000
CHEM_POT = -36.7264
TMPR = 300
L_MAX = 21
continue_simulations = True

xyz_fname = 'MD{:}_out.xyz'

# First, for CASSANDRA GCMC the governing CASSANDRA object should be created
css = cassandra.Cassandra()
# Read the some pre-prepared CASSANDRA .inp parameters file that can be modified then programatically
gcmcPropsRead = css.read_input('my_props.inp')

# Setting the gas system
co2 = system.read_lammps('co2.lmps')

# I Reading and arraging data from previous MC step
gas_mc = cassandra.McSystem(co2, chem_pot=CHEM_POT, max_ins=MAX_INS, is_rigid=True)
out_folder = 'results'

# Creating fixed polymer system
fixed_sst = system.read_lammps('pim.lmps')
fixed_sst.wrap()
l = 1

# This is used to continue simulations using the most fresh files taken from the SAME out_folder
if continue_simulations:
    # Figure out what is the last iteration with finished CASSANDRA
    da = cassandra.DataAnalyzer(path=out_folder, mc_fname_mask='*.before_md.lmps')
    l = int(re.search('\A\d+', os.path.split(da.file_names[-1])[1]).group())
    fixed_ref_sst = system.read_lammps(os.path.join(da.file_names[-1]))

    # Clearing previous LAMMPS log file if it exists
    prev_log = os.path.join(out_folder, str(l) + '.md.log')
    if os.path.exists(prev_log):
        os.remove(prev_log)

    # Updating simulation system with the data from the latest lmps file
    for p in fixed_sst.particles:
        p.x = fixed_ref_sst.particles[p.tag].x
        p.y = fixed_ref_sst.particles[p.tag].y
        p.z = fixed_ref_sst.particles[p.tag].z
    fixed_sst.dim = fixed_ref_sst.dim

gcmcPropsRead['Run_Name'] = str(l) + '.gcmc'
gcmc = cassandra.GCMC(gas_mc, fixed_sst, out_folder=out_folder, props_file=str(l) + '.gcmc_props.inp', **gcmcPropsRead)
css.add_gcmc(gcmc)

if continue_simulations:
    gcmc.upd_simulation()
else:
    css.run()

time_step = 0.5
lngth = 2e+4
nonrig_group_name = 'nonrig_bod'
rig_group_name = 'rig_bod'
while l < L_MAX:
    # >>> 2N: MD (LAMMPS) step:
    sim_sst = gcmc.tot_sst
    sim_sst.write_lammps(os.path.join(out_folder, str(l) + '.before_md.lmps'))
    sim = lmps.Simulation(sim_sst, log=os.path.join(out_folder, str(l) + '.md.log'), cutoff=14.0)

    # custom definitions for the neighbour list updates
    sim.add_custom('neighbor 1.0 nsq \nneigh_modify once no every 1 delay 0 check yes')

    # adding group definitions to separate rigid and non-rigid bodies
    grp_tmpl = 'group {:} id {:}'
    sim.add_custom(grp_tmpl.format(nonrig_group_name, gcmc.group_by_id('nonrigid')[0]))
    sim.add_custom(grp_tmpl.format(rig_group_name, gcmc.group_by_id('rigid')[0]))

    # create the description of the molecular dynamics simulation
    tmp_md = lmps.MolecularDynamics(ensemble='npt', timestep=time_step, length=int(lngth), thermo=2500, pressure=5,
                                    temp=TMPR, new_v=True, scale_v=True)

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
    sim.run(np=4)

    # Updating the size of the fixed system from the MD simulations for the next MC
    fixed_sst.dim = sim.system.dim

    # >>> 2N + 1: MC (CASSANDRA) step
    new_xyz_file = xyz_fname.format(str(l))

    # right after the MD the new .xyz file to read start configuration from is saved
    sim.system.write_xyz(os.path.join(out_folder, new_xyz_file))

    gcmcPropsRead['Run_Name'] = str(l + 1) + '.gcmc'

    # The 'Start_Type' property should be changed in order to make CASSANDRA read previously simulated system
    gcmcPropsRead['Start_Type'] = OrderedDict([('start_type', 'read_config'),
                                               ('species', [1] + [0] * len(gas_mc.chem_pot)),
                                               ('file_name', os.path.join(out_folder, new_xyz_file))])

    gas_mc = cassandra.McSystem(co2, chem_pot=CHEM_POT, max_ins=MAX_INS, is_rigid=True)
    gcmc = cassandra.GCMC(gas_mc, fixed_sst, out_folder=out_folder,
                          props_file=str(l + 1) + '.gcmc_props.inp', **gcmcPropsRead)
    css.add_gcmc(gcmc, is_new=True)
    css.run()
    l += 1
