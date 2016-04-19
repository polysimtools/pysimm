# ******************************************************************************
# pysimm.__init__ module
# ******************************************************************************
#
# error printing
# default executable path for lammps
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2016 Michael E. Fortunato
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

from __future__ import print_function
import os
from distutils.spawn import find_executable

__all__ = ['calc', 'forcefield', 'system', 'lmps', 'hoomd', 'utils', 'gasteiger', 'apps']

error = True
warning = True
verbose = True
debug = True

error_print = lambda *a, **k: print('(error) PySIMM:', *a) if error else lambda *a, **k: None
warning_print = lambda *a, **k: print('(warning) PySIMM:', *a) if warning else lambda *a, **k: None
verbose_print = lambda *a, **k: print('PySIMM:', *a) if verbose else lambda *a, **k: None
debug_print = lambda *a, **k: print('(debug) PySIMM:', *a) if debug else lambda *a, **k: None

if not os.environ.get('LAMMPS_EXEC'):
    os.environ['LAMMPS_EXEC'] = find_executable('lmp_serial') or find_executable('lmp_mpi')

if not os.environ.get('LAMMPS_EXEC'):
    if os.path.isfile(os.path.join(os.environ.get('HOME'), 'bin', 'lmp_mpi')):
        os.environ['LAMMPS_EXEC'] = os.path.join(os.environ.get('HOME'), 'bin', 'lmp_mpi')
    elif os.path.isfile('/usr/bin/lmp_mpi'):
        os.environ['LAMMPS_EXEC'] = '/usr/bin/lmp_mpi'
    else:
        os.environ['LAMMPS_EXEC'] = ('/apps/share64/debian7/lammps/'
                                     'lammps-06Apr15/bin/'
                                     'lmp_serial')
