# ******************************************************************************
# pysimm.cassandra module
# ******************************************************************************
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2017 Alexander Demidov, Michael E. Fortunato, Coray M. Colina
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from StringIO import StringIO
from subprocess import call, Popen, PIPE
import os
import re
import numpy as np
import random
import logging
import types
from collections import Iterable, OrderedDict
from pysimm import system
from string import ascii_uppercase
from pydoc import locate

DATA_PATH = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../dat/csndra_data'))

KCALMOL_2_K = 503.22271716452

CASSANDRA_EXEC = os.environ.get('CASSANDRA_EXEC')

# Creating a logger instance and send its output to console 'deafault'
logging.basicConfig(level=logging.INFO, datefmt='%H:%M:%S',
                    format='%(asctime)s [%(levelname)s]: %(message)s')

DEFAULT_PARAMS = {
    'Temperature_Info': 300,
    'Pressure_Info': 1,
    'Rcutoff_Low': 0.1
}


class MCSimulation(object):
    """pysimm.cassandra.MCSimulation

    Object containing the settings and the logic necessary to partially set-up an abstract Monte Carlo simulation
    to be submitted to the CASSANDRA software. The object also will include the simulation results once the simulations
    are finished.

    Attributes:
        mc_sst (:class:`~pysimm.cassandra.McSystem`) : describes all molecules to be inserted by CASSANDRA
        init_sst (:class:`~pysimm.system.System`) : describes the optional initial fixed molecular configuration for MC
            simulations (default: empty cubic box with 1 nm side length). If the particles in the system are not
            attributed with the flag `is_fixed` all of them are considered to be fixed, and will be marked with this
            flag, otherwise all particles with is_fixed=False will be removed.

    Keyword Args:
        out_folder (str) : the relative path of the simulation results (all .dat, .mcf, as well as .chk, ... files will
            go there). If the folder does not exist it will be created with 0755 permissions.
        props_file (str) : the name of the  .inp file.

    Note:
        Other keyword arguments that are accepted are the GCMC simulation settings. The keywords of the settings
        are the same as they are described in CASSANDRA specification but without # symbol.

        **For example**: the keyword argument `Run_Name='my_simulation'` will set `#Run_Name` setting in CASSANDRA
        input file to `my_simulation` value

    Parameters:
        props (dictionary) : include all simulation settings to be written to the CASSANDRA .inp file
        input (str) : text stream that will be written to the CASSANDRA .inp file
        tot_sst (:class:`~pysimm.system.System`) : object containing the results of CASSANDRA simulations
    """
    def __init__(self, mc_sst=None, init_sst=None, **kwargs):
        global DATA_PATH

        # Initializing CASSANDRA input stream, empty at the beginning
        self.input = ''

        # Initializing dictionary that contains records that directly will be sent to the .inp file
        self.props = OrderedDict()

        self.logger = logging.getLogger('MC Simulation')

        # Reading default properties of the GCMC simulations
        def_dat = Cassandra(system.System()).read_input(os.path.join(DATA_PATH, 'mc_default.inp'))

        tmp = kwargs.get('out_folder')  # Folder for the results and temporary files
        if tmp:
            self.out_folder = tmp
            if os.path.isabs(tmp):
                self.out_folder = os.path.relpath(tmp)
        else:
            self.out_folder = os.getcwd()
        if not os.path.exists(self.out_folder):
            os.makedirs(self.out_folder, mode=0755)
        prefix = kwargs.get('Run_Name', def_dat['Run_Name'])
        self.props['Run_Name'] = InpSpec('Run_Name', os.path.join(self.out_folder, prefix), '')

        self.props_file = os.path.join(self.out_folder, kwargs.get('props_file', ''))

        # Simple (one-value) dynamic properties
        self.props['Temperature_Info'] = InpSpec('Temperature_Info',
                                                 kwargs.get('Temperature_Info'), DEFAULT_PARAMS['Temperature_Info'])
        self.props['Pair_Energy'] = InpSpec('Pair_Energy', kwargs.get('Pair_Energy'), def_dat['Pair_Energy'])
        self.props['Rcutoff_Low'] = InpSpec('Rcutoff_Low', kwargs.get('Rcutoff_Low'), def_dat['Rcutoff_Low'])
        self.props['Mixing_Rule'] = InpSpec('Mixing_Rule', kwargs.get('Mixing_Rule'), def_dat['Mixing_Rule'])

        self.props['Seed_Info'] = InpSpec('Seed_Info', kwargs.get('Seed_Info'),
                                          [random.randint(int(1e+7), int(1e+8 - 1)),
                                           random.randint(int(1e+7), int(1e+8 - 1))])

        # Multiple-value one/many line dynamic properties
        self.props['Run_Type'] = InpSpec('Run_Type', kwargs.get('Run_Type'), def_dat['Run_Type'])
        self.props['Charge_Style'] = InpSpec('Charge_Style', kwargs.get('Charge_Style'), def_dat['Charge_Style'])
        self.props['VDW_Style'] = InpSpec('VDW_Style', kwargs.get('VDW_Style'), def_dat['VDW_Style'])
        self.props['Simulation_Length_Info'] = InpSpec('Simulation_Length_Info', kwargs.get('Simulation_Length_Info'),
                                                       def_dat['Simulation_Length_Info'],
                                                       **{'write_headers': True, 'new_line': True})
        self.props['CBMC_Info'] = InpSpec('CBMC_Info', kwargs.get('CBMC_Info'), def_dat['CBMC_Info'],
                                          **{'write_headers': True, 'new_line': True})

        self.props['Box_Info'] = InpSpec('Box_Info', kwargs.get('Box_Info'), def_dat['Box_Info'], **{'new_line': True})
        self.props['Property_Info 1'] = InpSpec('Property_Info 1', kwargs.get('Property_Info'), None, **{'new_line': True})

        # Setting the simulation total system
        if init_sst:
            self.tot_sst = init_sst.copy()
            self.tot_sst.center('box', [0, 0, 0], True)  # the center of the calculation box should be at origin
        else:
            self.logger.warning('The frame generating system for Monte-Carlo simulations is not set. '
                                'Creating empty cubic box of 1 nm size')
            self.tot_sst = system.System()
            self.tot_sst.forcefield = 'trappe/amber'
            self.tot_sst.dim = system.Dimension(dx=10, dy=10, dz=10)

        # Molecule configuration files describing all species of the system.
        # They are **absolutely** needed to start calculation
        mol_files = OrderedDict()

        # Some necessary verification of obtained system
        # TODO: check the forcefield to be sure that it is claas 1
        if False:
            self.logger.error('CASSANDRA supports only 1-st class force fields')
            exit(1)
        self.tot_sst.zero_charge()  # the sum of the charges should necessary be 0

        # Creating the system of fixed molecules
        self.fxd_sst_mcfile = None
        self.fxd_sst = kwargs.get('fixed_sst')
        if self.tot_sst.particles:
            tmp = self.tot_sst.copy()
            for p in tmp.particles:
                if not p.is_fixed:
                    tmp.particles.remove(p.tag)
            tmp.remove_spare_bonding()
            self.fxd_sst = tmp
            self.fxd_sst_mcfile = os.path.join(self.out_folder, 'fixed_syst.mcf')
            mol_files['file1'] = [self.fxd_sst_mcfile, 1]

        # Setting up the Monte Carlo system
        self.mc_sst = mc_sst
        if mc_sst:
            mc_sst.file_store = self.out_folder
            mol_files = mc_sst.update_props(mol_files)

        if kwargs.get('Molecule_Files'):
            mol_files = OrderedDict(sorted(kwargs.get('Molecule_Files').items()))

        # Raising an error and stop execution if no MCF information in one or another way is provided
        if (mc_sst is None) and (not kwargs.get('Molecule_Files')):
            self.logger.error('The molecular configuration files of gas molecules for simulation are not set. '
                              'Nothing to simulate. Exiting...')
            exit(0)

        self._n_spec = len(mol_files)
        self.props['Nbr_Species'] = InpSpec('Nbr_Species', self._n_spec, self._n_spec)
        self.props['Molecule_Files'] = InpSpec('Molecule_Files', mol_files, None, **{'new_line': True})

        # Synchronzing "start type" .inp record
        self.fxd_sst_xyz = ''
        pops_list = [0] * self._n_spec
        start_type = 'make_config'
        if self.fxd_sst:
            pops_list[0] = 1
            self.fxd_sst_xyz = os.path.join(self.out_folder, 'fixed_syst.xyz')
            start_type = 'read_config'
        start_conf_dict = OrderedDict([('start_type', start_type), ('species', pops_list),
                                       ('file_name', self.fxd_sst_xyz)])
        self.props['Start_Type'] = InpSpec('Start_Type', kwargs.get('Start_Type'), start_conf_dict)

        # Synchronzing Fragment files:
        frag_files = OrderedDict()
        if mc_sst:
            mc_sst.temperature = self.props['Temperature_Info'].value
            frag_files = mc_sst.update_frag_record(frag_files)
        if kwargs.get('Fragment_Files'):
            frag_files = OrderedDict(sorted(kwargs.get('Fragment_Files').items()))
        if (mc_sst is None) and (not kwargs.get('Fragment_Files')):
            self.logger.error('Cannot set the fragment files of gas molecules for simulation')
            exit(1)
        self.props['Fragment_Files'] = InpSpec('Fragment_Files', frag_files, None, **{'new_line': True})

    def write(self):
        """pysimm.cassandra.MCSimulation.write

        Iterates through the :class:`~MCSimulation.props` dictionary creating the text for correct CASSANDRA input
        """

        for key in self.props.keys():
            if self.props[key].value is not None:
                self.input += '{:}\n'.format(self.props[key].to_string())

        self.input += '\nEND'
        # Initializing output stream
        self.logger.info('Writing CASSANDRA .inp file to "{:}"...'.format(self.props_file))
        out_stream = open(self.props_file, 'w')
        out_stream.write('{:}'.format(self.input))
        out_stream.close()
        self.logger.info('File: "{:}" was created sucsessfully'.format(self.props_file))

    def group_by_id(self, group_key='matrix'):
        """pysimm.cassandra.MCSimulation.group_by_id

        Method groups the atoms of the system :class:`~MCSimulation.tot_sst` by a certain property. Will iterate through
        all atoms in the system and return indexes of only those atoms that match the property. Currently supports 3
        properties defined by the input keyword argument argument.

        Keyword Args:
            group_key (str): text constant defines the property to match. Possible keywords are:

                (1) `matrix` -- (default) indexes of the atoms in :obj:`~MCSimulation.fxd_sst`

                (2) `rigid` -- indexes of all atoms that have rigid atomic bonds. It is assumed here that rigid and
                    nonrigid atoms can interact only through intermolecular forces

                (3) `nonrigid` -- opposite of previous, indexes of all atoms that have nonrigid atomic bonds

        Returns:
            str:
                string in format `a1:b1 a2:b2 ...` where all indexes inside `[ak, bk]` belongs to the selected group
                and array of the form `[[a1, b1], [a2, b2], ...]`
        """
        fxd_sst_idxs = []
        if self.fxd_sst:
            fxd_sst_idxs = range(1, len(self.fxd_sst.particles) + 1)
        # Behaviour depending on type of particles to check
        check = lambda x: x
        if group_key.lower() == 'nonrigid':
            check = lambda x: not x.is_rigid
        elif group_key.lower() == 'rigid':
            check = lambda x: x.is_rigid
        elif group_key.lower() == 'matrix':
            check = lambda x: x.tag in fxd_sst_idxs
        idx_array = [[-1, -1]]
        for p in self.tot_sst.particles:
            if check(p):
                if idx_array[-1][0] > 0:
                    if abs(p.tag - idx_array[-1][1]) > 1:
                        idx_array.append([p.tag, p.tag])
                    else:
                        idx_array[-1][1] = p.tag
                else:
                    idx_array[-1] = [p.tag, p.tag]
        idx_string = ''
        for t in idx_array:
            if t[1] - t[0] > 1:
                idx_string += str(t[0]) + ':' + str(t[1]) + ' '
        return idx_string, idx_array

    def upd_simulation(self):
        """pysimm.cassandra.MCSimulation.upd_simulation

        Updates the :class:`~MCSimulation.tot_sst` field using the `MCSimulation.props['Run_Name'].chk` file. Will try
        to parse the checkpoint file and read the coordinates of the molecules inserted by CASSANDRA. If neither of the
        molecules from the :class:`~MCSimulation.mc_sst` can be fit to the text that was read the method will raise an
        exception. The fitting method: :class:`~McSystem.make_system` assumes that different molecules inserted by
        CASSANDRA have the same order of the atoms.
        """
        fname = '{:}{:}'.format(self.props['Run_Name'].value, '.chk')
        self.logger.info('Updating MC system from the CASSANDRA {:} file...'.format(fname))
        if os.path.isfile(fname):
            try:
                with open(fname, 'r') as inp:
                    lines = inp.read()
                    # Define the starting index of the lines with inserted atoms
                    start_ind = lines.find('total number of molecules')
                    end_ind = start_ind + lines[start_ind:-1].find('****', 1)
                    count_info = lines[start_ind:end_ind].split('\n')
                    offset = 1
                    if self.fxd_sst:
                        tmp = count_info[1].split()
                        offset += int(tmp[1]) * len(self.fxd_sst.particles)
                    # Grab the lines with inserted atoms
                    start_ind = lines.find('coordinates for all the boxes')
                    all_coord_lines = lines[start_ind:-1].split('\n')
                    inp.close()

                gas_lines = all_coord_lines[offset:]
                if len(gas_lines) > 0:
                    if self.fxd_sst:
                        self.tot_sst = self.fxd_sst.copy()
                    self.tot_sst.add(self.mc_sst.make_system(gas_lines), change_dim=False)
                    self.logger.info('Simulation system successfully updated')
                else:
                    self.logger.info('Final MC configuration has 0 new particles the initial system remains the same')

            except IndexError:
                self.logger.error('Cannot fit the molecules from the CASSANDRA file to the PySIMM system')
        else:
            self.logger.error('Cannot find the CASSANDRA checkpoint file to update simulation. '
                              'Probably it cannot be written by CASSANDRA to the place you specified')

    def __check_params__(self):
        """pysimm.cassandra.MCSimulation.__check_params__

        Private method designed for update the fields of the simulation object to make them conformed with each other
        """
        # Sync the simulation box parameters
        dx, dy, dz = self.tot_sst.dim.size()
        if (dx == dy) and (dy == dz):
            box_type = 'cubic'
            box_dims = str(dx)
        else:
            box_type = 'orthogonal'
            box_dims = '{0:} {1:} {2:}'.format(dx, dy, dz)

        upd_vals = OrderedDict([('box_count', 1),
                                ('box_type', box_type),
                                ('box_size', box_dims)])
        if ('Box_Info' in self.props.keys()) and isinstance(self.props['Box_Info'], InpSpec):
            self.props['Box_Info'] = InpSpec('Box_Info', upd_vals, None, **{'new_line': True})
        else:
            self.props['Box_Info'] = upd_vals

        tmp = self.props['Box_Info'].value['box_size'].split()
        if self.props['Box_Info'].value['box_type'] == 'cubic':
            tmp = tmp + tmp + tmp
        self.tot_sst.dim = system.Dimension(dx=float(tmp[0]), dy=float(tmp[1]), dz=float(tmp[2]))

        # Sync of the volume change frequency in equilibration regime
        if 'Prob_Volume' in self.props.keys():
            if self.props['Prob_Volume'] is None:
                self.props['Run_Type'].value['steps'] = self.props['Run_Type'].value['steps'][0]

    def __write_chk__(self, out_file):
        """pysimm.cassandra.MCSimulation.__write_chk__

        Creates the CASSANDRA checkpoint file basing on the information from the `~MCSimulation.tot_sst` field
        """
        # Initializing output stream
        if out_file == 'string':
            out_stream = StringIO()
        else:
            out_stream = open(out_file, 'w+')
        blk_separ = ' {:*^75}\n'

        # Writing Translation/rotation/... info
        out_stream.write(blk_separ.format('Translation,rotation, dihedral, angle distortion'))
        tmplate = '{t[0]$$}{t[1]$$}{t[2]$$}{t[3]$$}{t[4]$$}\n'
        molecules = self.props['Molecule_Files'].value
        for m, i in zip(molecules, range(len(molecules))):
            out_stream.write(tmplate.replace('$$', ':>6d').format(t=[i + 1, 0, 0, 0, 0]))
            out_stream.write(tmplate.replace('$$', ':>6d').format(t=[i + 1, 0, 0, 0, 0]))
            out_stream.write('{t[0]:>23.14E}{t[2]:>23.14E}{t[2]:>23.14E}\n'.format(t=[0, 0, 0]))
            out_stream.write('{0:>12d}{0:>12d}\n'.format(0, 0))

        # Small section with total # of MC trials -- it is 0 at the beginning
        out_stream.write(blk_separ.format('# of MC steps'))
        out_stream.write('{:>12d}\n'.format(0))

        # Writing Box-info information
        out_stream.write(blk_separ.format('Box info'))
        tmp = self.props['Box_Info'].value['box_size']
        x, y, z = 0, 0, 0
        bx_type = None
        if isinstance(tmp, types.ListType):
            if len(tmp) > 3:
                x, y, z = tmp[0], tmp[1], tmp[2]
        elif isinstance(tmp, int) or isinstance(tmp, float):
            x, y, z = tmp, tmp, tmp
        else:
            exit(0)

        # First 0 here correspond to the # of trials
        out_stream.write('{0:>12d}\n{1:<18.10f}\n{2:}\n'.format(0, x * y * z, self.props['Box_Info'].value['box_type']))

        tmpl = '{t[0]&&}{t[1]&&}{t[2]&&}\n'
        tmp = np.diag([x, y, z])
        for lines in tmp:
            out_stream.write((tmpl.replace('&&', ':^22.14f')).format(t=lines))

        tmp = np.diag([1 / x, 1 / y, 1 / z])
        for lines in tmp:
            out_stream.write((tmpl.replace('&&', ':^22.8f')).format(t=lines))
        out_stream.write('{:>18.12f}\n'.format(0))

        # Creating seeds
        out_stream.write(blk_separ.format('SEEDS'))
        out_stream.write('{t[0]:>12d}{t[1]:>12d}{t[2]:>12d}\n{t[3]:>12d}{t[4]:>12d}\n'.format(
            t=np.random.random_integers(int(1e+7), int(1e+8 - 1), 5)))

        # Writing total number of molecules by species
        out_stream.write(blk_separ.format('Info for total number of molecules'))
        out_stream.write('{0:>11d}{1:>11d}\n'.format(1, 1))  # Currentely only one polymer "molecule" in the simulation
        for i in range(1, len(molecules)):
            out_stream.write('{0:>11d}{1:>11d}\n'.format(i + 1, 0))

        out_stream.write(blk_separ.format('Writing coordinates of all boxes'))
        # Writing coordinates of atoms in all boxes
        line_template = '{l[0]:<5}{l[1]:<25.15f}{l[2]:<25.15f}{l[3]:<25.15f}{l[4]:>10d}\n'
        for parts in self.tot_sst.particles:
            try:
                out_stream.write(line_template.format(l=[parts.type.name, parts.x, parts.y, parts.z, 1]))
            except:
                continue
        out_stream.close()


class GCMC(MCSimulation):
    """pysimm.cassandra.GCMC
    Initiates the specific type of Monte Carlo simulations for CASSANDRA: simulations using Grand-Canonical ensemble of
    particles (constant volume-temperature-chemical potential, muVT). See :class:`~pysimm.cassandra.MCSimulation`
    for the detailed description of the properties.

    """
    def __init__(self, mc_sst=None, init_sst=None, **kwargs):
        MCSimulation.__init__(self, mc_sst, init_sst, **kwargs)

        self.logger.name = 'GCMC'
        self.props['Sim_Type'] = InpSpec('Sim_Type', 'GCMC', 'gcmc')

        # Path for all intermediate Cassandra files and results
        self.props_file = os.path.join(self.out_folder, kwargs.get('props_file', 'gcmc_input.inp'))

        add = 0
        if self.fxd_sst and self.fxd_sst.particles.count:
            add = 1
        self.props['Chemical_Potential_Info'] = InpSpec('Chemical_Potential_Info', kwargs.get('chem_pot'),
                                                        -30 * (self._n_spec - add))

        # Order of the next four items is IMPORTANT! Check the CASSANDRA spec file for further info
        def_init_prob = 0.25
        limits = [0.3] * self._n_spec
        if self.fxd_sst:
            limits[0] = 0
        self.props['Prob_Translation'] = InpProbSpec('Prob_Translation', kwargs.get('Prob_Translation'),
                                                     OrderedDict([('tot_prob', def_init_prob),
                                                                  ('limit_vals', limits)]),
                                                     **{'new_line': True, 'indicator': 'start'})
        tps = ['cbmc'] * self._n_spec
        if self.fxd_sst:
            tps[0] = 'none'
        self.props['Prob_Insertion'] = InpProbSpec('Prob_Insertion', kwargs.get('Prob_Insertion'),
                                                   OrderedDict([('tot_prob', def_init_prob), ('types', tps)]),
                                                   **{'new_line': True})

        self.props['Prob_Deletion'] = InpProbSpec('Prob_Deletion', kwargs.get('Prob_Deletion'), def_init_prob)

        max_ang = [180] * self._n_spec
        if self.fxd_sst:
            max_ang[0] = 0
        self.props['Prob_Rotation'] = InpProbSpec('Prob_Rotation', kwargs.get('Prob_Rotation'),
                                                  OrderedDict([('tot_prob', def_init_prob), ('limit_vals', max_ang)]),
                                                  **{'new_line': True,  'indicator': 'end'})


class NVT(MCSimulation):
    """pysimm.cassandra.NVT
    Initiates the specific type of Monte Carlo simulations for CASSANDRA: simulations using Canonical ensemble of
    particles (constant volume-temperature-number of particles, NVT). See :class:`~pysimm.cassandra.MCSimulation`
    for the detailed description of the properties.

    """
    def __init__(self, mc_sst=None, init_sst=None, **kwargs):
        MCSimulation.__init__(self, mc_sst, init_sst, **kwargs)
        self.logger.name = 'NVT'
        self.props_file = os.path.join(self.out_folder, kwargs.get('props_file', 'nvt-mc_input.inp'))
        self.props['Sim_Type'] = InpSpec('Sim_Type', 'nvt_mc', 'nvt_mc')

        move_probs = [1, 1, 1]
        limits = [0.3] * self._n_spec
        if self.fxd_sst:
            limits[0] = 0
        self.props['Prob_Translation'] = InpProbSpec('Prob_Translation', kwargs.get('Prob_Translation'),
                                                     OrderedDict([('tot_prob', move_probs[0]),
                                                                  ('limit_vals', limits)]),
                                                     **{'new_line': True, 'indicator': 'start'})
        sub_probs = [1] * self._n_spec
        if self.fxd_sst:
            sub_probs[0] = 0
        sm = sum(sub_probs)
        sub_probs = [s / sm for s in sub_probs]
        self.props['Prob_Regrowth'] = InpProbSpec('Prob_Regrowth', kwargs.get('Prob_Regrowth'),
                                                  OrderedDict([('tot_prob', move_probs[1]), ('sub_probs', sub_probs)]),
                                                  **{'new_line': True})
        max_ang = [180] * self._n_spec
        if self.fxd_sst:
            max_ang[0] = 0
        self.props['Prob_Rotation'] = InpProbSpec('Prob_Rotation', kwargs.get('Prob_Rotation'),
                                                  OrderedDict([('tot_prob', move_probs[2]), ('limit_vals', max_ang)]),
                                                  **{'new_line': True, 'indicator': 'end'})


class NPT(MCSimulation):
    """pysimm.cassandra.NPT
    Initiates the specific type of Monte Carlo simulations for CASSANDRA: simulations using Isobaric-Isothermal ensemble
    of particles (NPT). See :class:`~pysimm.cassandra.MCSimulation` for the detailed description of the properties.

    """
    def __init__(self, mc_sst=None, init_sst=None, **kwargs):
        MCSimulation.__init__(self, mc_sst, init_sst, **kwargs)

        # Initialising object attributes
        self.logger.name = 'NPT'
        self.props_file = os.path.join(self.out_folder, kwargs.get('props_file', 'npt-mc_input.inp'))

        # Initialising simulation-specific props attribute
        self.props['Sim_Type'] = InpSpec('Sim_Type', 'npt_mc', 'npt_mc')

        self.props['Pressure_Info'] = InpSpec('Pressure_Info',
                                              kwargs.get('Pressure_Info'), DEFAULT_PARAMS['Pressure_Info'])

        move_probs = [.34, .02, .32, .32]

        limits = [0.3] * self._n_spec
        if self.fxd_sst:
            limits[0] = 0
        self.props['Prob_Translation'] = InpProbSpec('Prob_Translation', kwargs.get('Prob_Translation'),
                                                     OrderedDict([('tot_prob', move_probs[0]),
                                                                  ('limit_vals', limits)]),
                                                     **{'new_line': True, 'indicator': 'start'})
        vol_margins = 0.1 * self.props['Box_Info'].value['box_size']
        self.props['Prob_Volume'] = InpProbSpec('Prob_Volume', kwargs.get('Prob_Volume'),
                                                OrderedDict([('tot_prob', move_probs[1]), ('types', vol_margins)]),
                                                **{'new_line': True})
        sub_probs = [1] * self._n_spec
        if self.fxd_sst:
            sub_probs[0] = 0
        sm = sum(sub_probs)
        sub_probs = [s / sm for s in sub_probs]
        self.props['Prob_Regrowth'] = InpProbSpec('Prob_Regrowth', kwargs.get('Prob_Regrowth'),
                                                  OrderedDict([('tot_prob', move_probs[2]), ('sub_probs', sub_probs)]),
                                                  **{'new_line': True})
        max_ang = [180] * self._n_spec
        if self.fxd_sst:
            max_ang[0] = 0
        self.props['Prob_Rotation'] = InpProbSpec('Prob_Rotation', kwargs.get('Prob_Rotation'),
                                                  OrderedDict([('tot_prob', move_probs[3]), ('limit_vals', max_ang)]),
                                                  **{'new_line': True, 'indicator': 'end'})


class InpSpec(object):
    """pysimm.cassandra.InpSpec

     Represents the most common object used for carrying one logical unit of the CASSANDRA simulation options

     Parameters:
         key (str) : the keyword of the simulation option (literally the string that goes after the # sign in
            CASSANDRA .inp file)
         value (object) : numerical or text values of the particular simulation option structured in a certain way.
            Here goes only the values that are wished to be changed (it might be just one field of a big dictionary)
         default (object) : the most complete default description of the simulation option

     Keyword Args:
         write_headers (boolean): if the :obj:`~value` is dictionary defines whether the dictionary keys should be
            written to the output
         new_line (boolean): if the :obj:`~value` is iterable defines whether each new element will be written to
            the new line
    """

    def __init__(self, key, value, default, **kwargs):
        self.key = key
        self.write_headers = kwargs.get('write_headers')
        self.is_new_line = kwargs.get('new_line')

        self.value = value
        if value:
            if isinstance(default, types.DictType):
                # Add from default structure all properties that were not defined by user
                for ky in value.keys():
                    default[ky] = value[ky]
                self.value = default
            else:
                self.value = value
        elif value == []:
            self.value = []
        else:
            # If nothing was passed write default
            self.value = default

    def to_string(self):
        """pysimm.cassandra.InpSpec.to_string

        Creates the proper text representation of the property stored in the :obj:`~value` field

        Returns:
            str:
                formatted text string
        """
        if self.value is not None:
            result = '# {:}\n'.format(self.key)
            # Strings
            if isinstance(self.value, types.StringTypes):
                result += str(self.value)
            # Dictionaries
            elif isinstance(self.value, types.DictType):
                for ks in list(self.value.keys()):
                    if self.write_headers:
                        result += ks + '  '

                    tmp = self.value[ks]
                    if (isinstance(tmp, Iterable)) & (not isinstance(tmp, types.StringTypes)):
                        result += '   '.join(str(p) for p in tmp)
                    else:
                        result += str(tmp)

                    if self.is_new_line:
                        result += '\n'
                    else:
                        result += ' '
                result = result[:-1]  # Remove the very last new line character
            # Lists
            elif isinstance(self.value, Iterable):
                for elem in self.value:
                    if isinstance(elem, Iterable):
                        subresult = ''
                        for subelem in elem:
                            subresult = subresult + str(subelem) + ' '
                    else:
                        subresult = str(elem) + ' '
                    result += subresult
            # Simple types
            else:
                result += str(self.value)
            result += '\n!{:^^20}\n'.format('')
            return result


class InpProbSpec(InpSpec):
    """pysimm.cassandra.InpSpec

    Extension of the :class:`~InpSpec` class that takes into account special representation of the movement
    probabilities in the CASSANDRA input file.
    """
    def __init__(self, key, value, default, **kwargs):
        super(InpProbSpec, self).__init__(key, value, default, **kwargs)

    def to_string(self):
        tmp = super(InpProbSpec, self).to_string()
        if self.key == 'Prob_Translation':
            tmp = '# Move_Probability_Info\n\n' + tmp
        elif self.key == 'Prob_Rotation':
            tmp += '\n# Done_Probability_Info\n'
        return tmp


class McSystem(object):
    """pysimm.cassandra.McSystem

    Wrapper around the list of :class:`~pysimm.system.System` objects. Each element in the list represents single
    molecule of a different specie that will be used during MC simulations. Additionally, the object is responsible for
    creating .dat and .mcf files needed for the simulation and reading back the CASSANDRA simulation results.

    Attributes:
        sst (list of :class:`~pysimm.system.System`) : items representing single molecules of different species to be
            inserted by CASSANDRA. If the sst is a list (not a single value) it is assumed that all of the following
            properties are synchronized with it by indexes.
        chem_pot (list of int) : chemical potential for each specie [Joule/mol]

    Keyword Args:
        max_ins (list of int) : defines the highest possible number of molecules of corresponding specie.
            Basing on these values CASSANDRA allocates memory for simulations. (default: 5000).
        is_rigid (list of boolean): defines whether the atoms in the particular molecule should be marked as rigid
            or not. **Important!** In current implementation the module doesn't support flexible molecule angles, so
            the `is_rigid=False` is designed to be used exclusively for **single bead** molecules.

    Parameters:
        made_ins (list of int) : number of particles of each specie inserted by CASSANDRA.
        mcf_file (list of str) : defines full relative names of molecule configuration files **(.mcf)** required by
            CASSANDRA. Files will be created automatically.
        frag_file (list of str) : defines full relative names of possible relative configuration files **(.dat)**
            required by CASSANDRA. Files will be created automatically.
    """
    def __init__(self, sst, **kwargs):
        self.logger = logging.getLogger('MC_SYSTEM')
        self.sst = make_iterable(sst)
        for sst in self.sst:
            # Checking that the force-field of the input system is of the class-1 as it is direct CASSANDRA restriction
            if isinstance(sst, system.System):
                sst.zero_charge()
                sst.add_particle_bonding()
                if sst.ff_class:
                    if not (sst.ff_class == '1'):
                        self.logger.error('Currently cassandra supports only with **Type-I** force fields. '
                                          'The PYSIMM systems you provided are of the different types'
                                          'Exiting...')
                        exit(1)
                else:
                    self.logger.info('The Force-Field type of the system is not defined. '
                                     'Assuming it is **Type-1** force field')

                    sst.ff_class = '1'

                if not all([pt.name for pt in sst.particle_types]):
                    self.logger.error('The name of at least one particle type in MC system is not defined. '
                                      'Will not be able to map particles back after the CASSANDRA simulations. '
                                      '\nPlease, setup the names for all particle types for your MC system')
                    exit(1)

            # Decorating the system with bonds_fixed flag and angle_fixed flag
            for bt in sst.bond_types:
                bt.is_fixed = True
            for at in sst.angle_types:
                if at.k > 70:
                    at.is_fixed = True

        self.file_store = os.getcwd()
        self.max_ins = make_iterable(kwargs.get('max_ins', 5000))
        self.is_rigid = make_iterable(kwargs.get('is_rigid', [True] * len(self.sst)))
        self.made_ins = [0] * len(self.sst)
        self.mcf_file = []
        self.frag_file = []
        self.temperature = None

    def update_props(self, props):
        """pysimm.cassandra.McSystem.update_props

        For each specie in the system creates the .mcf file required for CASSANDRA simulation.

        Args:
            props (dictionary) : contains the .mcf file names and maximally allowed number of molecules insertions.
                The dictionary is to be assigned to 'Molecule_Files' property of the MC simulation

        Returns:
            props: updated input dictionary
       """
        # Generate correct .mcf files
        al_ind = 0
        for (sstm, count) in zip(self.sst, range(len(self.sst))):
            fullfile = os.path.join(self.file_store, '{:}{:}{:}'.format('particle', str(count + 1), '.mcf'))
            for p_type in sstm.particle_types:
                if p_type.elem and (not p_type.real_elem):
                    p_type.real_elem = p_type.elem
                p_type.elem = ascii_uppercase[int(al_ind / 10)] + str(al_ind % 10)
                al_ind += 1
            McfWriter(sstm, fullfile).write()
            self.mcf_file.append(fullfile)
        # Make the files list to be returned
        offset = len(props)
        for (mcf, ins, count) in zip(self.mcf_file, self.max_ins, range(1 + offset, len(self.mcf_file) + 1 + offset)):
            props['file' + str(count)] = [mcf, ins]
        return props

    def update_frag_record(self, frag_record):
        """pysimm.cassandra.McSystem.update_frag_record

        For each specie in the system creates the single configuration .dat file required for CASSANDRA simulation.

        Args:
            frag_record: dictionary containing the .dat file names and their ids. The dictionary is to be assigned to
            'Molecule_Files' property of the MC simulation

        Returns:
            dictionary:
                updated dictionary
        """
        # Generating the structure files
        if self.temperature is None:
            self.temperature = 300

        for (sstm, count) in zip(self.sst, range(len(self.sst))):
            fullfile = os.path.join(self.file_store, '{:}{:}{:}'.format('particle', str(count + 1), '.dat'))
            with open(fullfile, 'w') as out:
                frag_count = 1
                out.write('{:>12d}\n'.format(frag_count))
                out.write('{:>21f}{:>21f}\n'.format(self.temperature, 0))
                tmplte = '{:<10}{:<24f}{:<24f}{:<24f}\n'
                for prt in sstm.particles:
                    out.write(tmplte.format(prt.type.elem, prt.x, prt.y, prt.z))
            self.frag_file.append(fullfile)
        # Generating the files list
        for (frags, count) in zip(self.frag_file, range(1, len(self.frag_file) + 1)):
            frag_record['file' + str(count)] = [frags, count]
        return frag_record

    def make_system(self, text_output):
        """pysimm.cassandra.McSystem.make_system

        Parses the checkpoint (.chk) file made by CASSANDRA and creates new molecules basing on the new coordinates
        information. Assumes that all atoms of a certain molecule are listed in .chk file together (molecule
        identifiers are not mixed).

        Note:
            The logic of comparison of the xyz-like text record from the .chk file with the
            :class:`~pysimm.system.System` object is most straightforward: It is the consecutive comparison of particle
            names and first letters (before the white space) in the text record. In this implementation order matters!
            For example, for CO2, if in the system atoms are ordered as C-O-O and in the text they are ordered as
            O-C-O fit will fail.

        Args:
            text_output (str): text stream from the CASSANDRA .chk file containing the coordinates of newly inserted
                molecules

        Returns:
            :class:`~pysimm.system.System` : object containing all newly inserted molecules
        """
        tmp_sst = None
        count = 0  # counter of the lines in the input file
        sys_idx = 0  # counter of the gas molecules to lookup
        while count < len(text_output):
            tmp = self.sst[sys_idx].copy()
            dictn = text_output[count:(len(tmp.particles) + count)]
            if self.__fit_atoms__(tmp, dictn):
                for p in tmp.particles:
                    vals = dictn[p.tag - 1].split()
                    # Read the coordinates from the text output of the CASSANDRA simulation
                    p.x, p.y, p.z = map(float, vals[1:4])
                    # Force velocities of the particles to be 0
                    p.vx, p.vy, p.vz = 0.0, 0.0, 0.0
                    p.molecule.syst_tag = 0
                if self.is_rigid[sys_idx]:
                    for p in tmp.particles:
                        p.is_rigid = True
                if tmp_sst:
                    tmp_sst.add(tmp)
                else:
                    tmp_sst = tmp.copy()
                self.made_ins[sys_idx] += 1
                count += len(tmp.particles)
                sys_idx = 0
            else:
                sys_idx += 1
                if sys_idx >= len(self.sst):
                    self.logger.error('Wasn\'t able to read CASSANDRA .chk file. '
                                      'Please check either MC-simulation provided to PySIMM or the CASSANDRA '
                                      'checkpoint file ')
                    exit(1)
        if tmp_sst:
            tmp_sst.update_tags()
            tmp_sst.objectify()
        return tmp_sst

    def __fit_atoms__(self, molec, text_lines):
        """pysimm.cassandra.McSystem.__fit_atoms__

        Implements simple logic of comparison of the xyz-like text record with the :class:`~pysimm.system.System`
        object. The comparison is based on the consecutive comparison of particle names and first letters (before the
        white space) in the text. In this implementation order matters!  E.g. for CO2, if in the system atoms are
        ordered as C-O-O and in the text they are ordered like O-C-O fit will return False.

        Returns:
            boolean: flag whether the text record fit the molecule or not
        """
        flag = True
        # Cannot map anything if number of molecules is different from number of data lines
        if len(molec.particles) != len(text_lines):
            return False
        # Check the sequence of element names they
        for p in molec.particles:
            vals = text_lines[p.tag - 1].split()

            if vals[0] != p.type.elem:
                return False
        return flag


class Cassandra(object):
    """pysimm.cassandra.Cassandra

    Organizational object for running CASSANDRA simulation tasks. In current implementation it is able to run Canonical,
    Grand Canonical, and Isothermal-Isobaric Monte Carlo simulations (:class:`~GCMC`, :class:`~NVT`, and :class:`~NPT`,
    correspondingly).

    Parameters:
        system (:class:`~pysimm.system.System`) : molecular updated during the simulations
        run_queue (list) : the list of scheduled tasks
    """
    def __init__(self, init_sst):
        self.logger = logging.getLogger('CSNDRA')

        # Assume all particles in initial system are fixed
        self.system = init_sst
        if init_sst.particles:
            for p in init_sst.particles:
                p.is_fixed = True
        self.run_queue = []

    def run(self):
        """pysimm.cassandra.Cassandra.run

        Method that triggers the simulations. Does two consecutive steps: **(1)** tries to write all files necessary
        for simulation (.dat, .inp, .mcf): **(2)** tries to invoke the CASSANDRA executable.

        """
        global CASSANDRA_EXEC

        if check_cs_exec():
            for task in self.run_queue:
                # Write .inp file
                task.write()
                # Write .xyz of the fixed system if provided
                if task.fxd_sst:
                    if task.fxd_sst_mcfile is not None:
                        McfWriter(task.fxd_sst, task.fxd_sst_mcfile).write('atoms')
                    task.fxd_sst.write_xyz(task.fxd_sst_xyz)
                try:
                    self.logger.info('Starting the GCMC simulations with CASSANDRA')
                    print('{:.^60}'.format(''))
                    p = Popen([CASSANDRA_EXEC, task.props_file], stdin=PIPE, stdout=PIPE, stderr=PIPE)
                    stout, sterr = p.communicate()
                    print(stout)
                    print(sterr)
                    task.upd_simulation()
                    self.system = task.tot_sst.copy()

                except OSError as ose:
                    self.logger.error('There was a problem calling CASSANDRA executable')
                    exit(1)
                except IOError as ioe:
                    if check_cs_exec():
                        self.logger.error('There was a problem running CASSANDRA. '
                                          'The process started but did not finish')
                        exit(1)
        else:
            self.logger.error('There was a problem running CASSANDRA: seems it is not configured properly.\n'
                              'Please, be sure the CSNDRA_EXEC environment variable is set to the correct '
                              'CASSANDRA executable path. The current path is set to:\n\n{}\n\n'.format(CASSANDRA_EXEC))
            exit(1)

    def add_simulation(self, ens_type, obj=None, **kwargs):
        """pysimm.cassandra.Cassandra.add_simulation

        Method for adding new Monte Carlo simulation to the run queue.

        Args:
            ens_type: Type of the molecular ensemble for the Monte-Carlo simulations. The supported options are: `GCMC`
                (Grand Canonical); `NVT` (canonical); `NPT` (isobaric-isothermal)
            obj: the entity that should be added. Will be ignored if it is not of a type  :class:`~MCSimulation`

        Keyword Args:
            is_new (boolean) : defines whether all previous simulations should be erased or not
            species (list of :class:`~pysimm.system.System`) : systems that describe molecules and will be passed to
                :class:`~McSystem` constructor.

        Note:
            Other keyword arguments of this method will be redirected to the :class:`~McSystem` and
            :class:`~MCSimulation` constructors. See their descriptions for the possible keyword options.
        """
        new_job = None

        # Reading the molecule ensemble type
        simul = locate('pysimm.cassandra.' + ens_type)
        if simul is None:
            self.logger.error('Unsopported simulation ensemble option. Please use ether GCMC, NPT, or '
                              'NVT in \'add_simulation\' ')
            exit(1)

        if isinstance(obj, MCSimulation):
            new_job = obj
        else:
            specs = kwargs.get('species')
            if specs:
                mc_sst = McSystem(specs, **kwargs)
                new_job = simul(mc_sst, self.system, **kwargs)
            else:
                self.logger.error('Incorrect ' + ens_type + ' initialization. Please provide either Cassandra.' +
                                  ens_type + ' simulation object or the dictionary with initialization parameters '
                                  'of that object')
                exit(1)
        # Clean the run queue if 'is_new' set to to True
        if kwargs.get('is_new'):
            self.run_queue[:] = []
        if new_job:
            new_job.__check_params__()
            self.run_queue.append(new_job)

    def add_gcmc(self, obj=None, **kwargs):
        """pysimm.cassandra.Cassandra.add_gcmc

        Ads new simulation in grand-canonical ensemble to the run queue.

        Args:
            obj: the entity that should be added. Will be ignored if it is not of a type  :class:`~GCMC`

        Keyword Args:
            is_new (boolean) : defines whether all previous simulations should be erased or not
            species (list of :class:`~pysimm.system.System`) : systems that describe molecules and will be passed to
                :class:`~McSystem` constructor.

        Note:
            Other keyword arguments of this method will be redirected to the :class:`~McSystem`, :class:`~MCSimulation`,
             and :class:`~GCMC` constructors. See their descriptions for the possible keyword options.
        """
        new_job = None
        if isinstance(obj, GCMC):
            new_job = obj
        else:
            specs = kwargs.get('species')
            if specs:
                mc_sst = McSystem(specs, **kwargs)
                new_job = GCMC(mc_sst, self.system, **kwargs)
            else:
                self.logger.error('Unknown GCMC initialization. Please provide either '
                                  'the dictionary with GCMC parameters or Cassandra.GCMC simulation object')
                exit(1)
        if kwargs.get('is_new'):
            self.run_queue[:] = []
        if new_job:
            new_job.__check_params__()
            self.run_queue.append(new_job)

    def add_npt_mc(self, obj=None, **kwargs):
        """pysimm.cassandra.Cassandra.add_npt_mc

        Ads new simulation in isobaric-isothermal ensemble to the run queue.

        Args:
            obj: the entity that should be added. Will be ignored if it is not of a type  :class:`~NPT`

        Keyword Args:
            is_new (boolean) : defines whether all previous simulations should be erased or not
            species (list of :class:`~pysimm.system.System`) : systems that describe molecules and will be passed to
                :class:`~McSystem` constructor.

        Note:
            Other keyword arguments of this method will be redirected to the :class:`~McSystem`, :class:`~MCSimulation`,
            and  :class:`~NPT` constructors. See their descriptions for the possible keyword options.
        """
        new_job = None
        if isinstance(obj, NPT):
            new_job = obj
        else:
            specs = kwargs.get('species')
            if specs:
                mc_sst = McSystem(specs, **kwargs)
                new_job = NPT(mc_sst, self.system, **kwargs)
            else:
                self.logger.error('Unknown NPT initialization. Please provide either '
                                  'the dictionary with NPT simulation parameters or Cassandra.NPT simulation object')
                exit(1)
        if kwargs.get('is_new'):
            self.run_queue[:] = []
        if new_job:
            new_job.__check_params__()
            self.run_queue.append(new_job)

    def add_nvt(self, obj=None, **kwargs):
        """pysimm.cassandra.Cassandra.add_nvt

        Ads new simulation in canonical ensemble to the run queue.

        Args:
            obj: the entity that should be added. Will be ignored if it is not of a type  :class:`~NVT`

        Keyword Args:
            is_new (boolean) : defines whether all previous simulations should be erased or not
            species (list of :class:`~pysimm.system.System`) : systems that describe molecules and will be passed to
                :class:`~McSystem` constructor.

        Note:
            Other keyword arguments of this method will be redirected to the :class:`~McSystem`, :class:`~MCSimulation`,
            and  :class:`~NVT` constructors. See their descriptions for the possible keyword options.
        """
        new_job = None
        if isinstance(obj, NVT):
            new_job = obj
        else:
            specs = kwargs.get('species')
            if specs:
                mc_sst = McSystem(specs, **kwargs)
                new_job = NVT(mc_sst, self.system, **kwargs)
            else:
                self.logger.error('Unknown NVT initialization. Please provide either '
                                  'the dictionary with NPT simulation parameters or Cassandra.NPT simulation object')
                exit(1)
        if kwargs.get('is_new'):
            self.run_queue[:] = []
        if new_job:
            new_job.__check_params__()
            self.run_queue.append(new_job)

    def read_input(self, inp_file):
        """pysimm.cassandra.Cassandra.read_input

        The method parses the CASSANDRA instructions file (.inp) split it into separate instructions and analyses each
        according to the instruction name.

        Args:
            inp_file (str) : the full relative path of the file to be read

        Returns:
            dictionary : read CASSANDRA properties in the format required by :class:`~GCMC`
        """
        result = {}
        if os.path.isfile(inp_file):
            self.logger.info('Reading simulation parameters from {:} file'.format(inp_file))
            # Reading the cassandra .inp file as one long string
            inp_stream = open(inp_file, 'r')
            lines = inp_stream.read()
            raw_props = lines.split('#')

            for prop in raw_props:
                line = re.sub('\n!.*', '', prop)  # Get rid of the CASSANDRA comments
                line = re.sub('\n(e|E)(n|N)(d|D)', '', line)  # Get rid of the 'END in the end of the file
                tmp = line.split()
                if len(tmp) > 1:
                    result[tmp[0]] = self.__parse_value__(tmp)

            # File seems fine let's close the stream and return true in the flag
            inp_stream.close()
            self.logger.info('Reading finished sucsessfully')
        else:
            self.logger.error('Cannot find specified file: \"{:}\"'.format(inp_file))
        return result

    def __parse_value__(self, cells):
        title = cells[0].lower()
        if title == 'run_type':
            return OrderedDict([('type', cells[1]), ('steps', map(int, cells[2:]))])

        elif title == 'charge_style':
            return OrderedDict([('type', cells[1]),
                                ('sum_type', cells[2]),
                                ('cut_val', float(cells[3])),
                                ('accuracy', float(cells[4]))])

        elif title == 'vdw_style':
            return OrderedDict([('type', cells[1]),
                                ('cut_type', cells[2]),
                                ('cut_val', float(cells[3]))])

        elif title == 'simulation_length_info':
            tmp = OrderedDict([('units', cells[2]),
                               ('prop_freq', int(cells[4])),
                               ('coord_freq', int(cells[6])),
                               ('run', int(cells[8]))])
            if len(cells) > 10:
                tmp['steps_per_sweep'] = int(cells[10])
                if len(cells) > 12:
                    tmp['block_averages'] = int(cells[12])
            return tmp

        elif title == 'cbmc_info':
            return OrderedDict([('kappa_ins', int(cells[2])),
                                ('kappa_dih', int(cells[4])),
                                ('rcut_cbmc', float(cells[6]))])

        elif title == 'box_info':
            size = float(cells[3])
            if len(cells) > 6:
                size = [float(cells[3]), float(cells[4]), float(cells[5])]
            return OrderedDict([('box_count', int(cells[1])), ('box_type', cells[2]), ('box_size', size)])

        elif title == 'prob_translation':
            vals = []
            for i in range(2, len(cells)):
                vals.append(float(cells[i]))
            return OrderedDict([('tot_prob', float(cells[1])),
                                ('limit_vals', vals)])

        elif title == 'prob_insertion':
            vals = []
            for i in range(2, len(cells)):
                vals.append(cells[i])
            return OrderedDict([('tot_prob', float(cells[1])),
                                ('types', vals)])

        elif title == 'prob_rotation':
            vals = []
            for i in range(2, len(cells)):
                vals.append(float(cells[i]))
            return OrderedDict([('tot_prob', float(cells[1])),
                                ('limit_vals', vals)])

        elif (title == 'molecule_files') or (title == 'fragment_files'):
            tmp = OrderedDict()
            for i, c in zip(range(1, len(cells) - 1, 2), range(1, 1 + len(cells) / 2)):
                tmp['file' + str(c)] = [cells[i], int(cells[i + 1])]
            return tmp

        elif title == 'start_type':
            if cells[1] == 'read_config':
                specs = []
                for i in range(2, len(cells) - 1):
                    specs.append(int(cells[i]))
                return OrderedDict([('start_type', 'read_config'),
                                    ('species', specs),
                                    ('file_name', cells[-1])])

            if cells[1] == 'make_config':
                specs = []
                for i in range(2, len(cells)):
                    specs.append(int(cells[i]))
                return OrderedDict([('start_type', 'make_config'),
                                    ('species', specs),
                                    ('file_name', '')])
            if cells[1] == 'add to config':
                self.logger.error('Sorry, \'add to config\' regime  of ''Start_Type option is not supported yet')
                exit(1)

            if cells[1] == 'checkpoint':
                self.logger.error('Sorry, \'checkpoint\' regime  of ''Start_Type option is not supported yet ')
                exit(1)

        elif title == 'property_info':
            if int(cells[1]) == 1:
                tmp = OrderedDict()
                for i in range(2, len(cells)):
                    tmp['prop' + str(i - 1)] = str.lower(cells[i])
                return tmp

        elif title == 'seed_info':
            return [int(cells[1]), int(cells[2])]

        elif (title == 'prob_deletion') or (title == 'rcutoff_low') or \
             (title == 'bond_prob_cutoff') or (title == 'chemical_potential_info'):
            return float(cells[1])

        elif (title == 'average_Info') or (title == 'nbr_species') or (title == 'temperature_info'):
            return int(cells[1])

        else:
            return cells[1]

    def unwrap_gas(self):
        """pysimm.cassandra.Cassandra.unwrap_gas

        Ensures that all particles that are not fixed are unwrapped, otherwise CASSANDRA might not interpret
        them correctly
        """
        gas_system = self.system.copy()
        for p in gas_system.particles:
            if p.is_fixed:
                gas_system.particles.remove(p.tag, update=False)
            else:
                self.system.particles.remove(p.tag, update=False)

        for m in gas_system.molecules:
            if any([t.is_fixed for t in m.particles]):
                gas_system.molecules.remove(m.tag, update=False)
            else:
                self.system.molecules.remove(m.tag, update=False)

        gas_system.remove_spare_bonding()
        self.system.remove_spare_bonding()
        gas_system.unwrap()
        self.system.add(gas_system, change_dim=False)


class McfWriter(object):
    """pysimm.cassandra.McfWriter

    Object responsible for creating the CASSANDRA Molecular Configuration file (.mcf).

    Attributes:
        syst (:class:`~pysimm.system.System`) :represents the molecule to be described
        file_ref (str) : full relative path to the file that will be created
    """
    # Section names in any .mcf file
    mcf_tags = ['# Atom_Info', '# Bond_Info', '# Angle_Info', '# Dihedral_Info',
                '# Improper_Info', '# Intra_Scaling', '# Fragment_Info', '# Fragment_Connectivity']
    empty_line = '0'

    def __init__(self, syst, file_ref):
        self.syst = syst
        self.file_ref = file_ref
        self.logger = logging.getLogger('MCF Writer')

    def write(self, typing='all'):
        """pysimm.cassandra.McfWriter.write

        Method creates the .mcf file writing only those sections of it that are marked to be written

        Args:
            typing (list) : the list of sections to be written or the text keyword. List items should be as they are
                defined in :class:`~pysimm.cassandra.McfWriter.mcf_tags` field); default 'all'
        """
        # Initializing output stream
        with open(self.file_ref, 'w') as out_stream:
            for (name, is_write) in zip(self.mcf_tags, self.__to_tags__(typing)):
                if is_write:
                    try:
                        method = getattr(self, '__write_' + str.lower(name[2:]) + '__')
                        method(out_stream)
                    except AttributeError:
                        self.__write_empty__(out_stream, name)
                else:
                    self.__write_empty__(out_stream, name)
            out_stream.write('\nEND')
            out_stream.close()

    def __write_empty__(self, out, name):
        out.write('{0:}\n{1:}\n\n'.format(name, self.empty_line))

    def __write_atom_info__(self, out):
        global KCALMOL_2_K
        text_tag = '# Atom_Info'
        if self.syst.particles.count > 0:
            # writing section header
            out.write('{:}\n'.format(text_tag))
            # Verify and fix net system charge
            self.syst.zero_charge()
            # writing total number of particles
            out.write('{0:<6}\n'.format(self.syst.particles.count))
            count = 0
            line_template = '{l[0]:<6}{l[1]:<7}{l[2]:<5}{l[3]:<8.3f}{l[4]:<10.6f}' \
                            '{l[5]:<6}{l[6]:<11.3f}{l[7]:<9.3f}\n'
            warn_flag = False
            for item in self.syst.particles:
                line = [count + 1, '', '', 0, 0, 'LJ', 0, 0]
                if item.charge:
                    line[4] = item.charge
                if item.type:
                    line[1] = item.type.tag
                    line[2] = item.type.tag
                    if item.type.name:
                        line[1] = item.type.name
                        line[2] = item.type.elem
                    else:
                        warn_flag = True
                    if item.type.mass:
                        line[3] = item.type.mass
                    if item.type.epsilon:
                        line[6] = KCALMOL_2_K * item.type.epsilon
                    if item.type.sigma:
                        line[7] = item.type.sigma
                else:
                    continue
                out.write(line_template.format(l=line))
                count += 1
            if warn_flag:
                self.logger.warning('Some particle type names (and/or element names) inside the system are not defined.'
                                    ' Will use type identifiers instead')
        else:
            self.__write_empty__(out, text_tag)
        out.write('\n')

    def __write_bond_info__(self, out):
        text_tag = '# Bond_Info'
        if self.syst.bonds.count > 0:
            # writing section header
            out.write('{:}\n'.format(text_tag))
            # writing total number of bonds
            out.write('{0:<6}\n'.format(self.syst.bonds.count))
            line_template = '{l[0]:<6d}{l[1]:<6d}{l[2]:<6d}{l[3]:<9}{l[4]:<6.3f}\n'
            count = 1
            for bond in self.syst.bonds:
                tmp = 'fixed'  # Fixed bond is the only option for CASSANDRA V-1.2
                line = [count, bond.a.tag, bond.b.tag, tmp, bond.type.r0]
                count += 1
                out.write(line_template.format(l=line))
            out.write('\n')
        else:
            self.__write_empty__(out, text_tag)

    def __write_angle_info__(self, out):
        text_tag = '# Angle_Info'
        if self.syst.angles.count > 0:
            # writing section header
            out.write('{:}\n'.format(text_tag))
            # writing total number of angles
            out.write('{0:<6}\n'.format(self.syst.angles.count))
            count = 1
            for angle in self.syst.angles:
                line_template = '{l[0]:<6d}{l[1]:<6d}{l[2]:<6d}{l[3]:<6d}{l[4]:<10}{l[5]:<13.3f}'
                line = [count, angle.a.tag, angle.b.tag, angle.c.tag]
                if hasattr(angle.type, 'is_fixed') and angle.type.is_fixed:
                    addon = ['fixed', angle.type.theta0]
                else:
                    addon = ['harmonic', KCALMOL_2_K * angle.type.k, angle.type.theta0]
                    line_template += '{l[6]:<13.3f}'
                count += 1
                out.write(line_template.format(l=line + addon) + '\n')
            out.write('\n')
        else:
            self.__write_empty__(out, text_tag)

    def __write_intra_scaling__(self, out):
        format_line = '{:<6.2f}{:<6.2f}{:<6.2f}{:<6.2f}'
        # writing section header
        out.write('{:}\n'.format('# Intra_Scaling'))
        # writing vdW scaling:  1-2 1-3 1-4 1-N
        out.write(format_line.format(0, 0, 0, 0) + '\n')
        # writing charge scaling:  1-2 1-3 1-4 1-N
        out.write(format_line.format(0, 0, 0, 0) + '\n\n')

    def __write_dihedral_info__(self, out):
        text_tag = '# Dihedral_Info'
        self.__write_empty__(out, text_tag)

    def __write_improper_info__(self, out):
        text_tag = '# Improper_Info'
        self.__write_empty__(out, text_tag)

    def __write_fragment_info__(self, out):
        # writing section header
        out.write('{:}\n'.format('# Fragment_Info'))
        # writing indexing
        out.write('{:}\n'.format(1))
        n = len(self.syst.particles)
        out.write('  '.join('{}'.format(item) for item in [1, n] + range(1, n + 1)))
        out.write('\n\n')

    def __write_fragment_connectivity__(self, out):
        text_tag = '# Fragment_Connectivity'
        self.__write_empty__(out, text_tag)

    def __to_tags__(self, inpt):
        n = len(self.mcf_tags)
        idxs = [True] * n
        if inpt.lower() == 'atoms':
            idxs = [False] * n
            idxs[self.mcf_tags.index('# Atom_Info')] = True
            idxs[self.mcf_tags.index('# Intra_Scaling')] = True
        return idxs


def check_cs_exec():
    """pysimm.cassandra.check_cs_exec

    Validates that the absolute path to the CASSANDRA executable is set in the `CASSANDRA_EXEC` environmental variable
    of the OS. The validation is called once inside the :class:`~Cassandra.run` method.
    """
    global CASSANDRA_EXEC
    flag = True
    if CASSANDRA_EXEC is None:
        print('Please specify the OS environment variable ''CASSANDRA_EXEC'' that points to '
              'CASSANDRA compiled binary file, which is by default cassandra_{compiler-name}[_openMP].exe ')
        flag = False
    return flag


def make_iterable(obj):
    """pysimm.cassandra.make_iterable

    Utility method that forces the attributes be iterable (wrap in a list if it contains of only one item)
    """
    it_obj = obj
    if not isinstance(obj, Iterable):
        it_obj = [obj]
    return it_obj

