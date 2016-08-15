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

ANTECHAMBER_EXEC  = os.environ.get('ANTECHAMBER_EXEC')

def calc_charges(s, charge_method='bcc'):
    ac = Antechamber(s)
    ac.convert_to_ac()
    ac.charges(charge_method)
    ac.cleanup()

class Antechamber(object):
    def __init__(self, s, **kwargs):
        self.system = s
        
    def convert_to_ac(self):
        self.system.write_pdb('pysimm.tmp.pdb', False)
        cl = '{} -fi pdb -i pysimm.tmp.pdb -fo ac -o pysimm.tmp.ac'.format(ANTECHAMBER_EXEC)
        p = Popen(cl.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
        p.communicate()
        
    def charges(self, charge_method='bcc'):
        cl = '{} -fi ac -i pysimm.tmp.ac -fo ac -o pysimm_charges.tmp.ac -c {}'.format(ANTECHAMBER_EXEC, charge_method)
        p = Popen(cl.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
        p.communicate()
        self.system.read_ac_charges('pysimm_charges.tmp.ac')
        
    def cleanup(self):
        fnames = ['pysimm.tmp.pdb', 'pysimm.tmp.ac', 'pysimm_charges.tmp.ac', ]
        fnames += ['ATOMTYPE.INF']
        fnames += glob.glob('ANTECHAMBER*')
        for fname in fnames:
            try:
                os.remove(fname)
            except:
                print('problem removing {}'.format(fname))
        