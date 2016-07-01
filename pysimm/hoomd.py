# ******************************************************************************
# pysimm.hoomd module
# ******************************************************************************
#
# in pre-development
# read_hoomd
# api to hoomd simulation code
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

import shlex
import shutil
import subprocess
import os
import re
import math
from math import sin,cos,pow,pi,floor,ceil
from random import randint
from random import gauss as gaussian
import numpy as np
from StringIO import StringIO
from xml.etree import ElementTree as ET

from pysimm import system
from pysimm import forcefield

from hoomd_script import *

nvt_integrator = None
npt_integrator = None
min_integrator = None


def read_hoomd(data_file, **kwargs):
    return system.read_hoomd(data_file, **kwargs)


def minimize(s, **kwargs):

    global hoomd_system, min_integrator

    name = kwargs.get('name') if kwargs.get('name') is not None else False
    log = kwargs.get('log') if kwargs.get('log') is not None else 'hoomd.log'
    reset = kwargs.get('reset') if kwargs.get('reset') is not None else False
    cutoff = kwargs.get('cutoff') if kwargs.get('cutoff') is not None else 12
    cutoff /= 10.
    timestep = kwargs.get('timestep') if kwargs.get('timestep') is not None else 1
    timestep /= 1000.
    ftol = kwargs.get('ftol') if kwargs.get('ftol') is not None else 1e-4
    etol = kwargs.get('etol') if kwargs.get('etol') is not None else 1e-4
    min_steps = kwargs.get('min_steps') if kwargs.get('min_steps') is not None else 10000
    write_dump = kwargs.get('dump')

    if reset:
        init.reset()
        hoomd_system = init_system(s, cutoff)
    elif not init.is_initialized():
        hoomd_system = init_system(s, cutoff)

    if min_integrator:
        min_integrator.enable()
        min_integrator.set_params(dt=timestep, ftol=ftol, Etol=etol)
    else:
        min_integrator = integrate.mode_minimize_fire(group=group.all(), dt=timestep, ftol=ftol, Etol=etol)

    if write_dump:
        dump.xml(filename='pst_hoomd_min.xml', vis=True)
        dump.dcd(filename='pst_hoomd_min.dcd', period=write_dump, overwrite=False, unwrap_full=True)
    run(min_steps)

    if min_integrator:
        min_integrator.disbale()

    if fire.has_converged():
        print 'converged'
    else:
        print 'not converged'

    for p in hoomd_system.particles:
        s.particles[p.tag+1].x = p.position[0] * 10
        s.particles[p.tag+1].y = p.position[1] * 10
        s.particles[p.tag+1].z = p.position[2] * 10
 
    if name:
        print '%s simulation using HOOMD successful' % name
    else:
        print 'molecular dynamics using HOOMD successful'


def md(s, **kwargs):

    global hoomd_system, nvt_integrator, npt_integrator

    name = kwargs.get('name') if kwargs.get('name') is not None else False
    log = kwargs.get('log') if kwargs.get('log') is not None else 'hoomd.log'
    reset = kwargs.get('reset') if kwargs.get('reset') is not None else False
    cutoff = kwargs.get('cutoff') if kwargs.get('cutoff') is not None else 12
    cutoff /= cutoff
    timestep = kwargs.get('timestep') if kwargs.get('timestep') is not None else 1
    timestep /= 1000
    ensemble = kwargs.get('ensemble') if kwargs.get('ensemble') is not None else 'nvt'
    temp = kwargs.get('temp') if kwargs.get('temp') is not None else 300
    temp *= 0.0083144621
    pressure = kwargs.get('pressure') if kwargs.get('pressure') is not None else 1
    new_v = kwargs.get('new_v') if kwargs.get('new_v') is not None else True
    length = kwargs.get('length') if kwargs.get('length') is not None else 2000
    write_dump = kwargs.get('dump')

    if reset:
        init.reset()
        hoomd_system = init_system(s, cutoff)
    elif not init.is_initialized():
        hoomd_system = init_system(s, cutoff)

    if new_v:
        px = py = pz = 0.0
        for p in hoomd_system.particles:
            mass = p.mass
            vx = gaussian(0, temp / mass)
            vy = gaussian(0, temp / mass)
            vz = gaussian(0, temp / mass)

            p.velocity = (vx, vy, vz)
            px += mass*vx
            py += mass*vy
            pz += mass*vz

        px /= s.nparticles
        py /= s.nparticles
        pz /= s.nparticles

        t = 0
        t_fixed = 0

        for p in hoomd_system.particles:
            mass = p.mass
            v = p.velocity
            t += (math.pow(v[0], 2)+math.pow(v[1], 2)+math.pow(v[2], 2))*p.mass
            p.velocity = (v[0] - px/mass, v[1] - py/mass, v[2] - pz/mass)
            t_fixed += (math.pow(v[0], 2)+math.pow(v[1], 2)+math.pow(v[2], 2))*p.mass

        print t/3/s.nparticles,t_fixed/3/s.nparticles

    all = group.all()

    integrate.mode_standard(dt=timestep)
    if ensemble == 'nvt':
        if nvt_integrator:
            nvt_integrator.enable()
            nvt_integrator.set_params(T=temp,tau=timestep*100)
        else:
            nvt_integrator = integrate.nvt(group=all,T=temp,tau=timestep*100)

    if write_dump:
        dump.xml(filename='pst_dump.xml',vis=True)
        dump.dcd(filename='pst_dump.dcd',period=write_dump,overwrite=True,unwrap_full=True)

    if log:
        analyze.log(filename=log,quantities=['pair_lj_energy', 'momentum','potential_energy','kinetic_energy','temperature'],period=100, header_prefix='#')

    run(length)

    if nvt_integrator: nvt_integrator.disable()
    if npt_integrator: npt_integrator.disable()


    for p in hoomd_system.particles:
        s.particles[p.tag+1].x = nm2a(p.position[0])
        s.particles[p.tag+1].y = nm2a(p.position[1])
        s.particles[p.tag+1].z = nm2a(p.position[2])

    if name: print '%s simulation using HOOMD successful'%name
    else: print 'molecular dynamics using HOOMD successful'


def init_system(s, nb_cut=1.2):

    s.write_hoomd('temp.xml')

    hoomd_system = init.read_xml(filename='temp.xml')

    os.remove('temp.xml')

    lj = pair.lj(r_cut=nb_cut)
    for pt1 in s.particle_types:
        for pt2 in s.particle_types:
            lj.pair_coeff.set(pt1.name, pt2.name, epsilon=math.sqrt(kcal2kj(pt1.epsilon)*kcal2kj(pt2.epsilon)),
                              sigma=a2nm(float(pt1.sigma)+float(pt2.sigma))/2)

    harmonic = bond.harmonic(name='harmonic_bonds')
    for bt in s.bond_types:
        harmonic.bond_coeff.set(bt.name,k=kcal2kj(bt.k)*2,r0=a2nm(bt.r0))

    harmonic_angle = angle.harmonic()
    for at in s.angle_types:
        harmonic_angle.set_coeff(at.name,k=kcal2kj(at.k)*2,t0=float(at.theta0)*math.pi/180)

    harmonic_dihedral = dihedral.harmonic()
    for dt in s.dihedral_types:
        harmonic_dihedral.set_coeff(dt.name, k=kcal2kj(dt.k)*2, d=int(dt.d), n=int(dt.n))

    return hoomd_system
