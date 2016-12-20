# ******************************************************************************
# pysimm.amber module
# ******************************************************************************
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2016 Michael E. Fortunato, Coray M. Colina
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

import os
import sys
import json
import glob
from subprocess import call, Popen, PIPE

from pysimm import forcefield
from pysimm import error_print
from pysimm import warning_print
from pysimm import debug_print


ANTECHAMBER_EXEC  = os.environ.get('ANTECHAMBER_EXEC')


def cleanup_antechamber():
    """pysimm.amber.cleanup_antechamber

    Removes temporary files created by antechamber and pysimm.

    Args:
        None

    Returns:
        None
    """
    fnames = ['pysimm.tmp.pdb', 'pysimm.tmp.ac' ]
    fnames += ['ATOMTYPE.INF']
    fnames += glob.glob('ANTECHAMBER*')
    for fname in fnames:
        try:
            os.remove(fname)
        except:
            print('problem removing {} during cleanup'.format(fname))


def calc_charges(s, charge_method='bcc', cleanup=True):
    """pysimm.amber.calc_charges

    Calculates charges using antechamber. Defaults to am1-bcc charges. 

    Args:
        s: System for which to calculate charges. System object is updated in place
        charge_method: name of charge derivation method to use (default: bcc)
        cleanup: removes temporary files created by antechamber (default: True)

    Returns:
        None
    """
    s.write_pdb('pysimm.tmp.pdb')
    cl = '{} -fi pdb -i pysimm.tmp.pdb -fo ac -o pysimm.tmp.ac -c {}'.format(ANTECHAMBER_EXEC, charge_method)
    p = Popen(cl.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.communicate()
    with file('pysimm.tmp.ac') as f:
        f.next()
        f.next()
        line = f.next()
        while line.split()[0] == 'ATOM':
            tag = int(line.split()[1])
            charge = float(line.split()[-2])
            s.particles[tag].charge = charge
            line = f.next()
    if cleanup:
        cleanup_antechamber()

        
def get_forcefield_types(s, types='gaff', f=None):
    """pysimm.amber.get_forcefield_types

    Uses antechamber to determine atom types. Defaults to GAFF atom types. Retrieves ParticleType objects from force field is provided 

    Args:
        s: System for which to type
        types: name of atom types to use (default: gaff)
        f: forcefield object to retrieve ParticleType objects from if not present in s (default: None)

    Returns:
        None
    """
    s.write_pdb('pysimm.tmp.pdb')
    cl = '{} -fi pdb -i pysimm.tmp.pdb -fo ac -o pysimm.tmp.ac -at {}'.format(ANTECHAMBER_EXEC, types)
    p = Popen(cl.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.communicate()
    with file('pysimm.tmp.ac') as fr:
        fr.next()
        fr.next()
        line = fr.next()
        while line.split()[0] == 'ATOM':
            tag = int(line.split()[1])
            type_name = line.split()[-1]
            if s.particle_types.get(type_name):
                s.particles[tag].type = s.particle_types.get(type_name)[0]
            elif f:
                pt = f.particle_types.get(type_name)
                if pt:
                    s.particles[tag].type = s.particle_types.add(pt[0].copy())
                    
            else:
                error_print('cannot find type {} in system or forcefield'.format(type_name))
            line = fr.next()
