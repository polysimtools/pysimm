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
        