from __future__ import print_function

import os
import sys
import imp
import shutil
from glob import glob
from pprint import pprint
from cStringIO import StringIO

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout
        
if os.path.isdir('tmp'):
    shutil.rmtree('tmp')
    os.mkdir('tmp')
else:
    os.mkdir('tmp')
os.chdir('tmp')

print('testing {:.<64}'.format('imports'), end='')
import pysimm
from pysimm import cli
from pysimm import amber
from pysimm import utils
from pysimm import calc
from pysimm import system
from pysimm import lmps
from pysimm import forcefield
from pysimm import apps
from pysimm import models

from pysimm.forcefield import *
from pysimm.apps import *
from pysimm.models.monomers.dreiding import *
from pysimm.models.monomers.gaff import *
from pysimm.models.monomers.gaff2 import *
print('passed')

example_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, 'Examples')

examples = glob(os.path.join(example_dir, '*'))

scripts = []
for example in examples:
    scripts += glob(os.path.join(example, '*.py'))
    scripts += glob(os.path.join(example, '*', '*.py'))

output = []
for script in sorted(scripts, key=lambda x: x.split('Examples/')[1]):
    print('testing Examples/{:.<55}'.format(script.split('Examples/')[1]), end='')
    foo = imp.load_source('create', script)
    try:
        with Capturing(output) as output:
            foo.run(test=True)
        print('passed')
    except:
        print('FAILED')
    
os.chdir('../')
shutil.rmtree('tmp')
