# ******************************************************************************
# pysimm.appps.poreblazer module
# ******************************************************************************
#
# api to poreblazer simulation code
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

import sys
from subprocess import Popen, call, PIPE
from random import randint
from time import strftime
import shlex

try:
    from Rappture.tools import getCommandOutput as RapptureExec
except ImportError:
    pass


def psd(s, **kwargs):
    """pysimm.apps.poreblazer.psd

    Perform pore size distribution calculation using PoreBlazer v2.0

    Args:
        atoms: file name to contain ff parameters (ff.atoms)
        data: file name to write xyz file (data.xyz)
        angles: angles of simlation box (90.0 90.0 90.0)
        insertions: number of insertions for calculation (500)
        min_probe: minimum probe size (1.0)
        probe_dr: step size to increase probe size (0.2)
        max_probe: maximum probe size: 25
        psd_save: T/F to save psd points (F)
        psd_range: range in which to save psd points (2.5,3.8)
        exec_path: path to poreblazer psd executable (psd.exe)
        gen_files: if True, only generate input do not execute (None)

    Returns:
        None
    """
    atoms = kwargs.get('atoms', 'ff.atoms')
    data = kwargs.get('data', "'data.xyz'")
    angles = kwargs.get('angles', '90.0 90.0 90.0')
    insertions = kwargs.get('insertions', 500)
    min_probe = kwargs.get('min_probe', 1.0)
    probe_dr = kwargs.get('probe_dr', 0.2)
    max_probe = kwargs.get('max_probe', 25)
    psd_save = kwargs.get('psd_save', 'F')
    psd_range = kwargs.get('psd_range', '2.5,3.8')

    exec_path = kwargs.get('exec_path', 'psd.exe')
    nanohub = kwargs.get('nanohub')

    gen_files = kwargs.get('gen_files')

    with open('psd.in', 'w+') as f:
        f.write('%s\n' % atoms)
        f.write('%s\n' % data)
        f.write('%s\n' % insertions)
        f.write('%s\n' % min_probe)
        f.write('%s\n' % probe_dr)
        f.write('%s\n' % max_probe)
        f.write('%s %s %s\n' % (s.dim.dx, s.dim.dy, s.dim.dz))
        f.write('%s\n' % angles)
        f.write('%s\n' % randint(10000, 99999))
        f.write('%s\n' % psd_save)
        f.write('%s\n' % psd_range)

    with open('ff.atoms', 'w+') as f:
        f.write('%s\n\n' % s.particle_types.count)
        for pt in s.particle_types:
            f.write('%s\t%f\n' % (pt.tag, pt.sigma))

    s.write_xyz(elem=False)

    if gen_files:
        return

    print('%s: starting pore size distribution simulation using poreblazer'
          % strftime('%H:%M:%S'))
    if nanohub:
        print('%s: sending pore size distribution simulation to computer cluster' % strftime('%H:%M:%S'))
        sys.stdout.flush()
        cmd = ('submit -n 1 -w %s -i psd.in -i ff.atoms -i data.xyz '
               'poreblazer-2.0.0_psd < psd.in'
               % (24*60))
        stdo, stde = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    else:
        stdin = open('psd.in')
        stdout = open('psd.out', 'w+')
        call(exec_path, stdin=stdin, stdout=stdout, shell=True)
        stdin.close()
        stdout.close()
    print('%s: pore size distribution simulation using poreblazer successful'
          % strftime('%H:%M:%S'))


def surface(s, **kwargs):
    """pysimm.apps.poreblazer.surface

    Perform accessible surface area calculation using PoreBlazer v2.0

    Args:
        atoms: file name to contain ff parameters (ff.atoms)
        data: file name to write xyz file (data.xyz)
        angles: angles of simlation box (90.0 90.0 90.0)
        insertions: number of insertions for calculation (1000)
        probe: probe size (3.681)
        probe_type: type of probe (hs)
        vis: True to save visual (F)
        exec_path: path to poreblazer surface executable (surface.exe)

    Returns:
        None
    """
    atoms = kwargs.get('atoms', 'ff.atoms')
    data = kwargs.get('data', "'data.xyz'")
    angles = kwargs.get('angles', '90.0 90.0 90.0')
    insertions = kwargs.get('insertions', 1000)
    probe = kwargs.get('probe', '3.681')
    probe_type = kwargs.get('probe_type', 'hs')
    vis = kwargs.get('vis') or 'F'

    exec_path = kwargs.get('exec_path', 'surface.exe')

    with open('surf_area.in', 'w+') as f:
        f.write('%s\n' % atoms)
        f.write('%s\n' % data)
        f.write('%s\n' % probe)
        f.write('%s\n' % insertions)
        f.write('%s %s %s\n' % (s.dim.dx, s.dim.dy, s.dim.dz))
        f.write('%s\n' % angles)
        f.write('%s\n' % probe_type)
        f.write('%s\n' % randint(10000, 99999))
        f.write('%s\n' % vis)

    with open('ff.atoms', 'w+') as f:
        f.write('%s\n\n' % s.particle_types.count)
        for pt in s.particle_types:
            f.write('%s\t%f\t%f\n' % (pt.tag, pt.sigma, pt.mass))

    s.write_xyz(elem=False)

    print('%s: starting surface area simulation using poreblazer'
          % strftime('%H:%M:%S'))
    stdin = open('surf_area.in')
    stdout = open('surf_area.out', 'w+')
    call(exec_path, stdin=stdin, stdout=stdout, shell=True)
    stdin.close()
    stdout.close()
    print('%s: surface area simulation using poreblazer successful'
          % strftime('%H:%M:%S'))

    s.surf_area = float(open('surf_area.out').readlines()[-1].split()[-1])

    return s.surf_area


def pore(s, **kwargs):
    """pysimm.apps.poreblazer.pore

    Perform pore volume calculation using PoreBlazer v2.0

    Args:
        atoms: file name to contain ff parameters (ff.atoms)
        data: file name to write xyz file (data.xyz)
        angles: angles of simlation box (90.0 90.0 90.0)
        insertions: number of insertions for calculation (1000)
        temp: temperature at which to perform simulation (300)
        pore_probe: sigma, epsilon, cutoff parameters for probe (2.58, 10.22, 12.8)
        exec_path: path to poreblazer pore executable (pore_he.exe)

    Returns:
        None
    """
    boltzmann_kcal = 0.001987204

    atoms = kwargs.get('atoms', 'ff.atoms')
    data = kwargs.get('data', "'data.xyz'")
    angles = kwargs.get('angles', '90.0 90.0 90.0')
    insertions = kwargs.get('insertions', 1000)
    temp = kwargs.get('temp', 300)
    pore_probe = kwargs.get('pore_probe', '2.58 10.22 12.8')

    exec_path = kwargs.get('exec_path', 'pore_he.exe')

    with open('pore_volume.in', 'w+') as f:
        f.write('%s\n' % atoms)
        f.write('%s\n' % data)
        f.write('%s\n' % insertions)
        f.write('%s %s %s\n' % (s.dim.dx, s.dim.dy, s.dim.dz))
        f.write('%s\n' % angles)
        f.write('%s\n' % temp)
        f.write('%s\n' % pore_probe)
        f.write('%s\n' % randint(10000, 99999))

    with open('ff.atoms', 'w+') as f:
        f.write('%s\n\n' % s.particle_types.count)
        for pt in s.particle_types:
            f.write('%s\t%f\t%f\t%f\n' % (pt.tag, pt.sigma,
                                          pt.epsilon/boltzmann_kcal, pt.mass))

    s.write_xyz(elem=False)

    print('%s: starting pore volume simulation using poreblazer'
          % strftime('%H:%M:%S'))
    stdin = open('pore_volume.in')
    stdout = open('pore_volume.out', 'w+')
    call(exec_path, stdin=stdin, stdout=stdout, shell=True)
    stdin.close()
    stdout.close()
    print('%s: pore volume simulation using poreblazer successful'
          % strftime('%H:%M:%S'))

    s.pore_volume = float(open('pore_volume.out').readlines()[-1].split()[-1])

    return s.pore_volume


def void(s, **kwargs):
    """pysimm.apps.poreblazer.void

    Perform pore volume calculation using PoreBlazer v2.0 assuming a probe size of 0 to calculate void volume

    Args:
        atoms: file name to contain ff parameters (ff.atoms)
        data: file name to write xyz file (data.xyz)
        angles: angles of simlation box (90.0 90.0 90.0)
        insertions: number of insertions for calculation (1000)
        temp: temperature at which to perform simulation (300)
        pore_probe: sigma, epsilon, cutoff parameters for probe (0.00, 10.22, 12.8)
        exec_path: path to poreblazer pore executable (pore_he.exe)

    Returns:
        None
    """
    boltzmann_kcal = 0.001987204

    atoms = kwargs.get('atoms', 'ff.atoms')
    data = kwargs.get('data', "'data.xyz'")
    angles = kwargs.get('angles', '90.0 90.0 90.0')
    insertions = kwargs.get('insertions', 1000)
    temp = kwargs.get('temp', 300)
    pore_probe = kwargs.get('pore_probe', '0.0 10.22 12.8')

    exec_path = kwargs.get('exec_path', 'pore_he.exe')

    with open('void_volume.in', 'w+') as f:
        f.write('%s\n' % atoms)
        f.write('%s\n' % data)
        f.write('%s\n' % insertions)
        f.write('%s %s %s\n' % (s.dim.dx, s.dim.dy, s.dim.dz))
        f.write('%s\n' % angles)
        f.write('%s\n' % temp)
        f.write('%s\n' % pore_probe)
        f.write('%s\n' % randint(10000, 99999))

    with open('ff.atoms', 'w+') as f:
        f.write('%s\n\n' % s.particle_types.count)
        for pt in s.particle_types:
            f.write('%s\t%f\t%f\t%f\n' % (pt.tag, pt.sigma,
                                          pt.epsilon/boltzmann_kcal, pt.mass))

    s.write_xyz(elem=False)

    print('%s: starting void volume simulation using poreblazer'
          % strftime('%H:%M:%S'))
    stdin = open('void_volume.in')
    stdout = open('void_volume.out', 'w+')
    call(exec_path, stdin=stdin, stdout=stdout, shell=True)
    stdin.close()
    stdout.close()
    print('%s: void volume simulation using poreblazer successful'
          % strftime('%H:%M:%S'))

    s.void_volume = float(open('void_volume.out').readlines()[-1].split()[-2])

    s.set_frac_free_volume()

    return s.void_volume, s.frac_free_volume
