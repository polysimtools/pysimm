from pysimm import system, lmps, cassandra
from collections import OrderedDict
import os
import re
import random


def mc_md(gas_sst, fixed_sst=None, **kwargs):
    """pysimm.apps.mc_md

    Performs the iterative hybrid Monte-Carlo/Molecular Dynamics (MC-MD) simulations using pysimm.lmps for MD and
    pysimm.cassandra for MD

    Args:
        gas_sst: list of pysimm.system.System objects each of which describes a different molecule to be inserted by MC
        fixed_sst: fixed during th MC steps group of atoms (default: None)
        mcmd_niter: number of MC-MD iteradions (default: 10)
        sim_folder: relative path to the folder with all simulation files (default: 'results')
        mc_props: dictionary describing all MC properties needed for simulations (see pysimm.cassandra.GCMC and
        pysimm.cassandra.GCMC.props for details)
        md_props: dictionary containing all Molecular Dynamics settings needed for simulations (see
        pysimm.lmps.Simulation and pysimm.lmps.MolecularDynamics for details)
    """

    nonrig_group_name = 'nonrigid_b'
    rig_group_name = 'rigid_b'
    n_iter = kwargs.get('mcmd_niter') or 10
    sim_folder = kwargs.get('sim_folder') or 'results'
    xyz_fname = os.path.join(sim_folder, 'MD{:}_out.xyz')
    l = 1

    # Creating fixed polymer system
    fs = None
    if fixed_sst:
        if isinstance(fixed_sst, str):
            fs = system.read_lammps(fixed_sst)
            fs.wrap()
        elif isinstance(fixed_sst, system.System):
            fs = fixed_sst
            fs.wrap()
        else:
            print('Cannot setup the fixed system for the simulations. Skipping this')

    # Set the one-molecule gas systems
    gases = []
    for g in cassandra.make_iterable(gas_sst):
        if isinstance(g, str):
            try:
                gases.append(system.read_lammps(g))
            except IOError:
                print('Cannot read file: {}\nExiting...'.format(g))
                exit(1)
        if isinstance(g, system.System):
            gases.append(g)

    if not gases:
        print('There are no gas molecules were specified correctely\nThe gas molecules are needed to start the '
              'MC-MD simulations\nExiting...')
        exit(1)

    css = cassandra.Cassandra(fixed_sst)
    # Set the Monte-Carlo properties:
    mcp = kwargs.get('mc_props')
    if mcp:
        CHEM_POT = cassandra.make_iterable(mcp.get('Chemical_Potential_Info'))
        if not CHEM_POT:
            print('Missing chemical potential info\nExiting...')
            exit(1)
    else:
        print('Missing the MC Simulation settings\nExiting...')
        exit(1)
    mcp['Start_Type'] = OrderedDict([('species', [1] + [0] * len(CHEM_POT))])

    # Set the Molecular-Dynamics properties:
    mdp = kwargs.get('md_props')
    if not mdp:
        print('Missing the MD Simulation settings\nExiting...')
        exit(1)

    while l < n_iter + 1:
        mcp['Run_Name'] = str(l) + '.gcmc'

        css.add_gcmc(species=gases, is_new=True, chem_pot=CHEM_POT,
                     is_rigid=mcp.get('rigid_type') or [False] * len(gases),
                     out_folder=sim_folder, props_file=str(l) + '.gcmc_props.inp', **mcp)
        css.run()

        # >>> 2N: MD (LAMMPS) step:
        sim_sst = css.final_sst
        sim_sst.write_lammps(os.path.join(sim_folder, str(l) + '.before_md.lmps'))
        sim = lmps.Simulation(sim_sst,
                              log=os.path.join(sim_folder, str(l) + '.md.log'),
                              print_to_screen=mdp.get('print_to_screen'),
                              cutoff=mdp.get('cutoff'))

        # custom definitions for the neighbour list updates
        sim.add_custom('neighbor 1.0 nsq \nneigh_modify once no every 1 delay 0 check yes')

        # adding group definitions to separate rigid and non-rigid bodies
        grp_tmpl = 'group {:} id {:}'
        sim.add_custom(grp_tmpl.format('matrix', css.run_queue[0].group_by_id('matrix')[0]))
        sim.add_custom(grp_tmpl.format(nonrig_group_name, css.run_queue[0].group_by_id('nonrigid')[0]))
        rigid_mols = css.run_queue[0].group_by_id('rigid')[0]
        if rigid_mols:
            sim.add_custom(grp_tmpl.format(rig_group_name, rigid_mols))

        # create the description of the molecular dynamics simulation
        tmp_md = lmps.MolecularDynamics(ensemble=mdp.get('ensemble'),
                                        timestep=mdp.get('timestep'),
                                        length=int(mdp.get('length')),
                                        thermo=mdp.get('thermo'),
                                        temp=mdp.get('temp'),
                                        pressure=mdp.get('pressure'),
                                        dump=int(mdp.get('dump')),
                                        dump_name=os.path.join(sim_folder, str(l) + '.md.dump'),
                                        scale_v=True)

        # obtain the simulation (LAMMPS) input in order to customly modify it later
        tmp_md.write(sim)

        # replace single default fix with two separate fixes for rigid and nonrigid bodies separately
        old_line = re.search('(?<=(\nfix)).*', tmp_md.input).group(0)
        corr_fix = re.sub('all', nonrig_group_name, old_line)

        if rigid_mols:
            corr_fix += ' dilate all\n'
        else:
            corr_fix += '\n'

        if rigid_mols:
            corr_fix += 'fix' + re.sub('iso\s+\d+[.\d]*\s+\d+[.\d]*\s+\d+[.\d]*', '', old_line).\
                        replace('1', '2', 1). \
                        replace('all', rig_group_name). \
                        replace('npt', 'rigid/nvt/small molecule') + '\n'

        # adding the spring fix to the geometrical center of the system to avoid system creep
        corr_fix += 'fix {:} {:} spring tether {:} {:} {:} {:} {:}\n'.format(3, 'matrix', 30.0, 0.0, 0.0, 0.0, 0.0)

        # saving all fixes to the input
        tmp_md.input = tmp_md.input.replace(old_line, corr_fix)

        # adding "run 0" line for correct temperature scaling of the system with rigid molecules
        tmp_md.input = tmp_md.input.replace('velocity all scale',
                                            'velocity all create {:} {:}\nrun 0\nvelocity all scale'
                                            .format(mdp.get('temp'), random.randint(int(1e+5), int(1e+6) - 1)))

        # The input for correct simulations is set, starting LAMMPS:
        sim.add_custom(tmp_md.input)
        sim.run(np=1)

        # Updating the size of the fixed system from the MD simulations and saving the coordinates for the next MC
        css.init_sst.dim = sim.system.dim
        sim.system.write_xyz(xyz_fname.format(l))
        mcp['Start_Type']['file_name'] = xyz_fname.format(l)
        mcp['Start_Type']['species'] = [1] + [0] * len(CHEM_POT)
        l += 1

    return sim.system