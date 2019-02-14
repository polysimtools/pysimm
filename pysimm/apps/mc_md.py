from pysimm import system, lmps, cassandra
from collections import OrderedDict
import os
import re
import random
import types


def mc_md(gas_sst, fixed_sst=None, mcmd_niter=None, sim_folder=None, mc_props=None, md_props=None, **kwargs):
    """pysimm.apps.mc_md

    Performs the iterative hybrid Monte-Carlo/Molecular Dynamics (MC/MD) simulations using :class:`~pysimm.lmps` for
    MD and :class:`~pysimm.cassandra` for MC

    Args:
        gas_sst (list of :class:`~pysimm.system.System`) : list items describe a different molecule to be
            inserted by MC
        fixed_sst (:class:`~pysimm.system.System`) : fixed during th MC steps group of atoms (default: None)


    Keyword Args:
        mcmd_niter (int) : number of MC-MD iterations (default: 10)
        sim_folder (str): relative path to the folder with all simulation files (default: 'results')
        mc_props (dictionary) : description of  all MC properties needed for simulations (see
            :class:`~pysimm.cassandra.GCMC` and :class:`~pysimm.cassandra.GCMC.props` for details)
        md_props (dictionary):  description of all Molecular Dynamics settings needed for simulations (see
            :class:`~pysimm.lmps.Simulation` and :class:`~pysimm.lmps.MolecularDynamics` for details)

    Returns:
        :class:`~pysimm.system.System`:
            Final state of the simulated system
    """

    nonrig_group_name = 'nonrigid_b'
    rig_group_name = 'rigid_b'
    n_iter = mcmd_niter or 10
    sim_folder = sim_folder or 'results'
    xyz_fname = os.path.join(sim_folder, '{:}.md_out.xyz')
    l = 1

    # Creating fixed polymer system
    fs = None
    if fixed_sst:
        if isinstance(fixed_sst, system.System):
            fs = fixed_sst
            fs.wrap()
        else:
            print('Cannot setup the fixed system for the simulations. Skipping this')

    # Set the one-molecule gas systems
    gases = []
    if gas_sst:
        if isinstance(gas_sst, system.System):
            gases = [gas_sst]
        elif isinstance(gas_sst, types.ListType):
            for g in cassandra.make_iterable(gas_sst):
                if isinstance(g, system.System):
                    gases.append(g)

    if not gases:
        print('There are no gas molecules were specified correctely\nThe gas molecules are needed to start the '
              'MC-MD simulations\nExiting...')
        exit(1)

    css = cassandra.Cassandra(fixed_sst)
    # Set the Monte-Carlo properties:
    mcp = mc_props
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
    sim = None
    mdp = md_props
    if not mdp:
        print('Missing the MD Simulation settings\nExiting...')
        exit(1)

    while l < n_iter + 1:
        # >>> MC (CASSANDRA) step:
        mcp['Run_Name'] = str(l) + '.gcmc'

        css.add_gcmc(species=gases, is_new=True, chem_pot=CHEM_POT,
                     is_rigid=mcp.get('rigid_type') or [False] * len(gases),
                     out_folder=sim_folder, props_file=str(l) + '.gcmc_props.inp', **mcp)
        css.run()

        # >>> MD (LAMMPS) step:
        sim_sst = css.system.copy()
        sim_sst.write_lammps(os.path.join(sim_folder, str(l) + '.before_md.lmps'))
        sim = lmps.Simulation(sim_sst, print_to_screen=mcp.get('print_to_screen', False),
                              log=os.path.join(sim_folder, str(l) + '.md.log'))

        sim.add(lmps.Init(cutoff=mdp.get('cutoff'),
                          special_bonds=mdp.get('special_bonds'),
                          pair_modify=mdp.get('pair_modify')))

        # custom definitions for the neighbour list updates
        sim.add_custom('neighbor 1.0 nsq \nneigh_modify once no every 1 delay 0 check yes')

        # adding group definitions to separate rigid and non-rigid bodies
        sim.add(lmps.Group('matrix', 'id', css.run_queue[0].group_by_id('matrix')[0]))
        sim.add(lmps.Group(nonrig_group_name, 'id', css.run_queue[0].group_by_id('nonrigid')[0]))
        rigid_mols = css.run_queue[0].group_by_id('rigid')[0]
        if rigid_mols:
            sim.add(lmps.Group(rig_group_name, 'id', rigid_mols))

        # adding "run 0" line before velocities rescale for correct temperature init of the system with rigid molecules
        sim.add(lmps.Velocity(style='create'))
        if rigid_mols:
            sim.add_custom('run 0')
            sim.add(lmps.Velocity(style='scale'))

        # create the description of the molecular dynamics simulation
        sim.add_md(lmps.MolecularDynamics(name='main_fix',
                                          group=nonrig_group_name if rigid_mols else 'all',
                                          ensemble='npt',
                                          timestep=mdp.get('timestep'),
                                          temperature=mdp.get('temp'),
                                          pressure=mdp.get('pressure'),
                                          run=False,
                                          extra_keywords={'dilate': 'all'} if rigid_mols else {}))

        # create the second NVT fix for rigid molecules that cannot be put in NPT fix
        if rigid_mols:
            sim.add(lmps.MolecularDynamics(name='rig_fix',
                                           ensemble='rigid/nvt/small molecule',
                                           timestep=mdp.get('timestep'),
                                           length=mdp.get('length'),
                                           group=rig_group_name,
                                           temperature=mdp.get('temp'),
                                           pressure=mdp.get('pressure'),
                                           run=False))

        # add the "spring tether" fix to the geometrical center of the system to avoid system creep
        sim.add_custom('fix tether_fix matrix spring tether 30.0 0.0 0.0 0.0 0.0')
        sim.add(lmps.OutputSettings(thermo=mdp.get('thermo'),
                                    dump={'filename': os.path.join(sim_folder, str(l) + '.md.dump'),
                                          'freq': int(mdp.get('dump'))}))
        sim.add_custom('run {:}\n'.format(mdp.get('length')))

        # The input for correct simulations is set, starting LAMMPS:
        sim.run(np=mdp.get('np', 1))

        # Updating the size of the fixed system from the MD simulations and saving the coordinates for the next MC
        # css.system.dim = sim.system.dim
        css.system = sim.system.copy()
        css.unwrap_gas()
        css.system.write_xyz(xyz_fname.format(l))

        mcp['Start_Type']['file_name'] = xyz_fname.format(l)
        mcp['Start_Type']['species'] = [1] + css.run_queue[-1].mc_sst.made_ins
        l += 1

    return sim.system if sim else None

