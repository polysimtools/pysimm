# ******************************************************************************
# pysimm.cli module
# ******************************************************************************
#
# command line interface for pysimm
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

from __future__ import print_function

import argparse

supported_forcefields = ['dreiding', 'pcff', 'gaff']

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Command line tools '
                                                 'for pysimm')

    parser.add_argument('-v', dest='verbosity',
                        help='verbosity level for output from pysimm tools\n'
                             '1: error output\n'
                             '2: error, warning output\n'
                             '3: error, warning, verbose output (default)\n'
                             '4: error, warning, verbose, debugging output',
                        type=int, default=3, choices=[1, 2, 3, 4])

    parser.add_argument('--unwrap', action='store_true',
                        help="unwrap system so bonds do not cross "
                             "periodic boundary")

    parser.add_argument('--lammps', dest='lammps_data',
                        help="read a lammps data file and "
                             "store in local variable 's'")

    parser.add_argument('--mol', dest='molfile',
                        help="read a mol file and "
                             "store in local variable 's'")

    parser.add_argument('--cml', dest='cml_file',
                        help="read a cml file and "
                             "store in local variable 's'")

    parser.add_argument('--yaml', dest='yaml_file',
                        help="read a yaml file and "
                             "store in local variable 's'")
                             
    parser.add_argument('--smiles', dest='smiles',
                        help="create system from smiles using pubchem API and "
                             "store in local variable 's'")

    parser.add_argument('--forcefield', dest='forcefield',
                        help='forcefield name to type system with')

    parser.add_argument('--lmps2xyz', dest='lmps2xyz',
                        nargs=2, help='convert lammps data file to xyz file')

    parser.add_argument('--lmps2cdjson', dest='lmps2cdjson',
                        nargs=2, help='convert lammps data file to chemdoodle json file')

    parser.add_argument('--vmd', action='store_true',
                        help="visualize using mvd")

    parser.add_argument('--pymol', action='store_true',
                        help="visualize using mvd")

    args = parser.parse_args()

    import pysimm
    if args.verbosity == 4:
        pysimm.debug = pysimm.verbose = pysimm.warning = pysimm.error = True
    elif args.verbosity == 3:
        pysimm.verbose = pysimm.warning = pysimm.error = True
    elif args.verbosity == 2:
        pysimm.warning = pysimm.error = True
    elif args.verbosity == 1:
        pysimm.error = True
    pysimm.error_print = lambda *a, **k: print('(error) PySIMM:', *a) if pysimm.error else lambda *a, **k: None
    pysimm.warning_print = lambda *a, **k: print('(warning) PySIMM:', *a) if pysimm.warning else lambda *a, **k: None
    pysimm.verbose_print = lambda *a, **k: print('PySIMM:', *a) if pysimm.verbose else lambda *a, **k: None
    pysimm.debug_print = lambda *a, **k: print('(debug) PySIMM:', *a) if pysimm.debug else lambda *a, **k: None

    print('\nWelcome to the pySIMM command line interface\n')
    print('This is no more than a python2.7 interactive shell with certain pySIMM modules imported for your convenience')
    print('Importing modules now...')

    from pysimm import system, amber, lmps, forcefield

    if args.lammps_data:
        s = system.read_lammps(args.lammps_data)
    elif args.molfile:
        s = system.read_mol(args.molfile)
        if args.forcefield and args.forcefield in supported_forcefields:
            if args.forcefield.lower() == 'dreiding':
                print('typing with %s' % args.forcefield)
                s.apply_forcefield(forcefield.Dreiding())
            elif args.forcefield.lower() == 'pcff':
                s.apply_forcefield(forcefield.Pcff())
        elif args.forcefield:
            print('forcefield %s is not supported in '
                  'command line interface at this time')
    elif args.cml_file:
        s = system.read_cml(args.cml_file)
    elif args.yaml_file:
        s = system.read_yaml(args.yaml_file)
    elif args.smiles:
        s = system.read_pubchem_smiles(args.smiles)
    elif args.lmps2xyz:
        s = system.read_lammps(args.lmps2xyz[0])
        if args.unwrap:
            print('unwrapping system so bonds do no cross simulation boundaries...'
                  'this may take a while if your system is large')
            s.unwrap()
        s.write_xyz(args.lmps2xyz[1])
    elif args.lmps2cdjson:
        s = system.read_lammps(args.lmps2cdjson[0])
        s.write_chemdoodle_json(args.lmps2cdjson[1])

    if args.unwrap:
        try:
            print('unwrapping system so bonds do no cross simulation boundaries...'
                  'this may take a while if your system is large')
            s.unwrap()
        except Exception as e:
            print('\nsomething went wrong')
            print(e.message)

    if args.vmd:
        s.viz()
    elif args.pymol:
        s.viz('pymol')

    print('\nHave fun with your molecular systems!\n')
