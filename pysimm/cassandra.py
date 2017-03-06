# ******************************************************************************
# pysimm.cassandra module
# ******************************************************************************
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)

from StringIO import StringIO
from pysimm import lmps
from pysimm.utils import Item, ItemContainer
import re
import os
import numpy as np
import logging
import types
from collections import Iterable, OrderedDict

class GCMC(object):

    def __init__(self, **kwargs):
        # Text output stream
        self.input = ''

        # Dictionary containing properties for all .inp lines
        self.props = OrderedDict()

        # Static (unchangable) properties
        self.props['Sim_Type'] = InpSpec('gcmc', 'gcmc')

        # Simple (one-value) dynamic properties
        self.props['Run_Name'] = InpSpec(kwargs.get('Run_Name'), 'result.out')
        self.props['Temperature_Info'] = InpSpec(kwargs.get('Temperature_Info'), 300)
        self.props['Chemical_Potential_Info'] = InpSpec(kwargs.get('Chemical_Potential_Info'), -10)
        self.props['Average_Info'] = InpSpec(kwargs.get('Average_Info'), 1)
        self.props['Pair_Energy'] = InpSpec(kwargs.get('Pair_Energy'), 'true')
        self.props['Rcutoff_Low'] = InpSpec(kwargs.get('Rcutoff_Low'), 0.0)
        self.props['Seed_Info'] = InpSpec(kwargs.get('Seed_Info'),
                                          np.random.random_integers(1e+7, 1e+8 - 1, [1, 2]))
        self.props['Mixing_Rule'] = InpSpec(kwargs.get('Mixing_Rule'), 'lb')
        self.props['Bond_Prob_Cutoff'] = InpSpec(kwargs.get('Bond_Prob_Cutoff'), 1e-10)


        # Multiple-value one/many line dynamic properties
        self.props['Run_Type'] = InpSpec(kwargs.get('Run_Type'),
                                         OrderedDict([('type','Equilibration'),
                                                      ('steps', 100)]))

        self.props['Charge_Style'] = InpSpec(kwargs.get('Charge_Style'),
                                             OrderedDict([('type','coul'),
                                                          ('sum_type', 'ewald'),
                                                          ('cut_val', 40.00),
                                                          ('accuracy', 1e-5)]))

        self.props['VDW_Style'] = InpSpec(kwargs.get('VDW_Style'),
                                          OrderedDict([('type', 'lj'),
                                                       ('cut_type', 'cut_tail'),
                                                       ('cut_val', 40.00)]))

        self.props['Simulation_Length_Info'] = InpSpec(kwargs.get('Simulation_Length_Info'),
                                                       OrderedDict([('units', 'steps'),
                                                                    ('prop_freq', 100),
                                                                    ('coord_freq', 1000),
                                                                    ('run', 10000)]),
                                                       **{'write_headers': True, 'new_line': True})
        self.props['CBMC_Info'] = InpSpec(kwargs.get('CBMC_Info'),
                                          OrderedDict([('kappa_ins', 12),
                                                       ('kappa_rot', 0),
                                                       ('kappa_dih', 12),
                                                       ('rcut_cbmc', 6.5)]),
                                          **{'write_headers': True, 'new_line': True})

        self.props['Box_Info'] = InpSpec(kwargs.get('Box_Info'),
                                         OrderedDict([('box_count', 1),
                                                      ('box_type', 'cubic'),
                                                      ('box_size', 100)]),
                                         **{'new_line': True})

        #   Here you might define order of properties written to the file once the order is important
        self.props['Prob_Translation'] = None
        self.props['Prob_Insertion'] = None
        self.props['Prob_Deletion'] = None

        # Get number of species from the description of .mcf files that user must provide in order to run simulations
        nSpec = len(kwargs.get('Molecule_Files'))
        self.props['Nbr_Species'] = InpSpec(kwargs.get('Nbr_Species'), nSpec)

        # Files needed to start calculation
        self.props['Molecule_Files'] = InpSpec(OrderedDict(sorted(kwargs.get('Molecule_Files').items())),
                                               None, **{'is_file': True, 'new_line': True})

        self.props['Fragment_Files'] = InpSpec(OrderedDict(sorted(kwargs.get('Fragment_Files').items())),
                                               None, **{'is_file': True, 'new_line': True})

        self.props['Start_Type'] = InpSpec(kwargs.get('Start_Type'),
                                           OrderedDict([('start_type', 'read_config'),
                                                        ('species', 1),
                                                        ('file_name', 'some_file.xyz')]),
                                               **{'is_file': True})

        self.props['Property_Info 1'] = InpSpec(kwargs.get('Property_Info'),
                                              None, **{'new_line': True})

        # Friking exception!!!
        self.props['Prob_Translation'] = InpSpec(kwargs.get('Prob_Translation'),
                                                 OrderedDict([('tot_prob', 0.4),
                                                              ('limit_vals', [0.0000, 0.3600])]),
                                                 **{'new_line': True})

        self.props['Prob_Insertion'] = InpSpec(kwargs.get('Prob_Insertion'),
                                                 OrderedDict([('tot_prob', 0.3),
                                                              ('types', ['none', 'cbmc'])]),
                                                 **{'new_line': True})

        self.props['Prob_Deletion'] = InpSpec(kwargs.get('Prob_Deletion'), 0.3)



    def write(self, out_file):
        for key in self.props.keys():
            if self.props[key].value is not None:

                # Crutch no. 1:
                if key == 'Prob_Translation':
                    self.input = self.input + '# Move_Probability_Info\n\n'

                self.input = self.input + '# {:}\n'.format(key)
                self.input = self.input + '{:}\n'.format(self.props[key].to_string())

                # Crutch no. 2:
                if key == 'Prob_Deletion':
                    self.input = self.input + '\n# Done_Probability_Info\n'

                self.input = self.input + '!{:-^20}\n\n'.format('')

        self.input = self.input + '\nEND'

        #TODO: move that wired stuff from here once all infrastructure is ready
        # Initializing output stream
        if out_file == 'string':
            out_stream = StringIO()
        else:
            out_stream = open(out_file, 'w+')
        out_stream.write('{:}'.format(self.input))



class InpSpec(object):
    def __init__(self, value, template, **kwargs):

        self.write_headers = kwargs.get('write_headers')
        self.is_new_line = kwargs.get('new_line')
        self.is_file = kwargs.get('is_file')

        if value:
            if isinstance(template, types.DictType):
                #Add from default structure all properties that were not defined by user
                for key in value.keys():
                    template[key] = value[key]
                self.value = template
            else:
                self.value = value
        else:
            # If nothing was passed write default
            self.value = template



    def to_string(self):
        result = ''
        #Strings
        if isinstance(self.value, types.StringTypes):
            result = str(self.value)
        #Dictionaries
        elif isinstance(self.value, types.DictType):
            for ks in list(self.value.keys()):
                if self.write_headers:
                    result = result + ks + '   '

                tmp = self.value[ks]
                if (isinstance(tmp, Iterable)) & (not isinstance(tmp, types.StringTypes)):
                    result = result + '   '.join(str(p) for p in tmp)
                else:
                    result = result + str(tmp)

                if self.is_new_line:
                    result = result + '\n'
                else:
                    result = result + ' '
            result = result[:-1] # Remove the very last new line character
        #Lists
        elif isinstance(self.value, Iterable):
            for elem in self.value:
                if isinstance(elem, Iterable):
                    subresult = ''
                    for subelem in elem:
                        subresult = subresult + str(subelem) + ' '
                else:
                    subresult = str(elem)
                result = subresult + '\n'
        #Simple types
        else:
            result = str(self.value)
        return result



class InpFileSpec(InpSpec):
    def __init__(self, value, template, **kwargs):
        super(InpFileSpec, self).__init__(value, template, **kwargs)

        # Check the existence of the file-type variables. Continue only when the files exist!
        names = template['file_name']
        if not isinstance(names, types.StringTypes):
            if isinstance(names, Iterable):
                for theName in names:

                    InpSpec()

                    if not os.path.isfile(theName):
                        print("ERROR: cannot find a file " + theName + ".\n Please specify the file.\n" + " Aborting execution")
                        sys.exit(0)


    def to_string(self):
        return ''



class SimBox(Item):
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)





class Cassandra(object):
    """
    pysimm.cassandra.Cassandra
    Organizational object for Cassandra simulation that is able to run
    e.g. Gibbs Canonical Monte-Carlo (GCMC) simulations (see the GCMC class)

    """

    def __init__(self, s, **kwargs):
        
        self.kcalMol2K = 503.22271716452
        self.system = s
        self.sim = []

        self.boxes = kwargs.get('boxes') or ItemContainer()
        
        self.__defineStatics__()
        
        #Creating a logger instance and send its output to console 'deafault'
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger('CSNDRA');


    def run(self):
        try:
            call_cs_exec(self, np, nanohub)
        except OSError as ose:
            raise PysimmError('There was a problem calling CASSANDRA executable'), None, sys.exc_info()[2]
        except IOError as ioe:
            if check_lmps_exec():
                raise PysimmError('There was a problem running CASSANDRA. The process started but did not finish '
                                  'successfully. Check the generated log file'), None, sys.exc_info()[2]
            else:
                raise PysimmError('There was a problem running LAMMPS. LAMMPS is not configured properly. '
                                  'Make sure the LAMMPS_EXEC environment variable is set to the correct LAMMPS '
                                  'executable path. The current path is set to:\n\n{}'.format(LAMMPS_EXEC)), None, sys.exc_info()[2]


    def add_gcmc(self, template=None, **kwargs):
        if template is None:
            self.sim.append(GCMC(**kwargs))
        elif isinstance(template, GCMC):
            self.sim.append(template)


    def write_mcf(self, out_file):
        
        #Default value to fill in unimportant .mcf file categories
        emptyVal = 0
        
        # Initializing output stream
        if out_file == 'string':
            out_stream = StringIO()
        else:
            out_stream = open(out_file, 'w+')

        # Writing "Static content"
        for tag in self.mcfTags:
            out_stream.write('{0:}\n{1:d}\n\n'.format(tag, emptyVal))

        # Writing important for .mcf atomic data
        sys = self.system # alias of the System object
        out_stream.write('{:}\n'.format(self.atm_info_str))

        # writing total number of particles
        out_stream.write('{0:5}\n'.format(sys.particles.count))

        # writing atom-specific data
        line_template = '{l[0]:>5}{l[1]:>7}{l[2]:>5}{l[3]:>3.0f}{l[4]:>14.9f}{l[5]:>3}{l[6]:>11.6f}{l[7]:>11.6f}\n'
        count = 0

        if sys.particles.count > 0:
            for parts in sys.particles:
                line = [count + 1,  '',  '',  '',  0, 'LJ', 0, 0]
                if hasattr(parts,  'charge'):
                    line[4] = parts.charge
                if hasattr(parts,  'type'):
                    if hasattr(parts.type,  'name'):
                        line[1] = parts.type.name
                    if hasattr(parts.type,  'elem'):
                        line[2] = parts.type.elem
                    if hasattr(parts.type,  'mass'):
                        line[3] = parts.type.mass
                    if hasattr(parts.type,  'epsilon'):
                        line[6] = self.kcalMol2K * parts.type.epsilon
                        line[7] = parts.type.sigma
                    else:
                        continue
                    out_stream.write(line_template.format( l=line ))
                else:
                    continue
                count = count + 1
            out_stream.write('\nEND')
        out_stream.close()


    def parseBoxes(self, cells):

        for i in range(1, len(cells), 2):
            if cells[i] == 'CUBIC':
                txt = cells[i + 1].split()
                x = float(txt[0])
                y = float(txt[(len(txt) - 1)/2])
                z = float(txt[len(txt) - 1])
                vol = x * y * z

                tmpDict = {'bxType': 'CUBIC', 'x': x, 'y': y, 'z': z, 'vol': vol}
                self.boxes.add(SimBox(**tmpDict))



    def write_chk(self, out_file):
        # Initializing output stream
        if out_file == 'string':
            out_stream = StringIO()
        else:
            out_stream = open(out_file, 'w+')
        sys = self.system # alias of the System object

        blkSepar = '{:*^75}\n'

        #Writing Translation/rotation/... info
        contNfo = self.props['# Molecule_Files']
        out_stream.write(blkSepar.format('Translation,rotation, dihedral, angle distortion'))
        tmplate = '{t[0]$$}{t[1]$$}{t[2]$$}{t[3]$$}{t[4]$$}\n'

        for i in range(len(contNfo)):
            out_stream.write(tmplate.replace('$$', ':>6d').format(t =  map(int, np.insert(np.zeros(4), 0, i + 1))))
            out_stream.write(tmplate.replace('$$', ':>6d').format(t =  map(int, np.insert(np.zeros(4), 0, i + 1))))
            #TODO There are some nonzeros in Tylangas .chk file for index 2; check where they come from
            out_stream.write('{t[0]:>23.14E}{t[2]:>23.14E}{t[2]:>23.14E}\n'.format(t =  np.zeros(3)))
            out_stream.write('{0:>12d}{0:>12d}\n'.format(0, 0))

        #Small section with total # of MC trials -- it is 0 at the beggining
        out_stream.write(blkSepar.format('# of MC steps'))
        out_stream.write('{:>12d}\n'.format(0))

        #Writing Box-info stuff
        out_stream.write(blkSepar.format('Box info'))
        for box in self.boxes:
            #First 0 in input correspond to the # of trials
            out_stream.write('{0:>12d}\n{1:<18.10f}\n{2:}\n'.format(0, box.vol, box.bxType))
            
            tmpl = '{t[0]&&}{t[1]&&}{t[2]&&}\n'
            tmp = np.diag( [box.x, box.y, box.z] )
            for lines in tmp:
                out_stream.write((tmpl.replace('&&', ':^22.14f')).format(t=lines))

            tmp = np.diag( [1/box.x, 1/box.y, 1/box.z])
            for lines in tmp:
                out_stream.write((tmpl.replace('&&', ':^22.8f')).format(t = lines))
            out_stream.write('{:>18.12f}\n'.format(0))#TODO: Maximal volume displacement


        #Writing SEEDS !!!111
        out_stream.write(blkSepar.format('SEEDS'))
        out_stream.write('{t[0]:>12d}{t[1]:>12d}{t[2]:>12d}\n{t[3]:>12d}{t[4]:>12d}\n'.format(
                        t = np.random.random_integers(1e+7, 1e+8 - 1, 5)))

        #Writing total number of molecules by species
        out_stream.write(blkSepar.format('Info for total number of molecules'))
        out_stream.write('{0:>11d}{1:>11d}\n'.format(1,1))#Currentely one polymer molecule in the simulation
        for i in range(1, len(contNfo)):
            out_stream.write('{0:>11d}{1:>11d}\n'.format(i + 1,0))

        out_stream.write(blkSepar.format('Writing coordinates of all boxes'))
        # Writing coordinates of atoms in all boxes
        lineTemplate = '{l[0]:>6} {l[1]:>13.8f} {l[2]:>13.8f} {l[3]:>13.8f}\n {l[4]:>12d} \n'
        for parts in sys.particles:
            line = ['',  0,  0,  0, 1] #TODO: change the "1" to the actual box identifier
            try:
                line[0] = parts.type.name
                line[1] = parts.x
                line[2] = parts.y
                line[3] = parts.z
            except:
                continue
            out_stream.write(lineTemplate.format( l=line ))
        out_stream.close()


    def readParams(self, paramsFile):
        result = 0
        if os.path.isfile(paramsFile):
            #File seems fine let's return true in the flag
            self.logger.info('Reading parameters from {:} file'.format(paramsFile))
            result = 1
            
            #Reading the cassandra .inp file as one long string
            inp_stream = open(paramsFile, 'r')
            lines = inp_stream.read()
            
            #Splitting the long string to parts by hash-tags
            lists = []
            for i in range(len(self.inpTags)):
                lists.append(self.__parseParamString__(lines, self.inpTags[i]))

            self.props = dict(zip(self.inpTags, lists))
            
            self.parseBoxes(self.props['# Box_Info'])
            
            inp_stream.close()
            self.logger.info('Reading finished sucsessfully')

        return result



    def __parseParamString__(self, strng, header):
        valsList = None
        headerStart = '#'
        carrReturn = '\n'
        
        tmpInd = strng.find(header)
        if tmpInd > 0:
            startInd = strng.find(carrReturn, tmpInd) + 1
            endInd = strng.find(headerStart, startInd) - 1
            tmp = strng[startInd:endInd].split('\n')
        return [i for i in tmp if i] #Deleting the empty strings



    def __defineStatics__(self):
        #Sections of the Cassandra .mcf output file
        self.mcfTags = ['# Bond_Info',
                        '# Angle_Info',
                        '# Dihedral_Info',
                        '# Fragment_Info',
                        '# Improper_Info',
                        '# Fragment_Connectivity']

        #Sections of the Cassandra .chk output file
        self.chkTags = ['Translation,rotation, dihedral, angle distortion',
                        'of MC steps',
                        'Box info',
                        'SEEDS',
                        'Info for total number of molecules',
                        'Writing coordinates for all the boxes']

        #Tags of the Cassandra .inp file that important in forming of the 
        # correct .chk and .mcf files
        self.inpTags = ['# Box_Info', 
                        '# Molecule_Files']

        self.atm_info_str = '# Atom_Info'













