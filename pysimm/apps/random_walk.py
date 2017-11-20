# ******************************************************************************
# pysimm.apps.random_walk module
# ******************************************************************************
#
# psuedo random walk algorithm written using pysimm tools
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


from time import strftime
from itertools import permutations, izip

import numpy as np

from pysimm import system, lmps, forcefield, calc
from pysimm import error_print


def find_last_backbone_vector(s, m):
    """pysimm.apps.random_walk.find_last_backbone_vector

    Finds vector between backbone atoms in terminal monomer. Requires current system s, and reference monomer m.

    Args:
        s: :class:`~pysimm.system.System` object
        m: :class:`~pysimm.system.System` object
    Returns:
        list of vector components
    """
    head_pos = [0, 0, 0]
    tail_pos = [0, 0, 0]
    for p in s.particles[-1*m.particles.count:]:
        if p.linker == 'head':
            head_pos = [p.x, p.y, p.z]
        elif p.linker == 'tail':
            tail_pos = [p.x, p.y, p.z]
    return [head_pos[0] - tail_pos[0], head_pos[1] - tail_pos[1], head_pos[2] - tail_pos[2]]


def copolymer(m, nmon, s_=None, **kwargs):
    """pysimm.apps.random_walk.copolymer

    Builds copolymer using random walk methodology using pattern

    Args:
        m: list of reference monomer :class:`~pysimm.system.System`s
        nmon: total number of monomers to add to chain
        s_: :class:`~pysimm.system.System` in which to build polymer chain (None)
        settings: dictionary of simulation settings
        density: density at which to build polymer (0.3)
        forcefield: :class:`~pysimm.forcefield.Forcefield` object to acquire new force field parameters
        capped: True/False if monomers are capped
        unwrap: True to unwrap final system
        traj: True to build xyz trajectory of polymer growth (True)
        pattern: list of pattern for monomer repeat units, should match length of m ([1 for _ in range(len(m))])
        limit: during MD, limit atomic displacement by this max value (LAMMPS ONLY)
        sim: :class:`~pysimm.lmps.Simulation` object for relaxation between polymer growth
    Returns:
        new copolymer :class:`~pysimm.system.System`
    """
    m = [x.copy() for x in m]

    settings = kwargs.get('settings', {})
    density = kwargs.get('density', 0.3)
    f = kwargs.get('forcefield')
    capped = kwargs.get('capped')
    unwrap = kwargs.get('unwrap')
    traj = kwargs.get('traj', True)
    pattern = kwargs.get('pattern', [1 for _ in range(len(m))])
    limit = kwargs.get('limit', 0.1)
    sim = kwargs.get('sim')

    for m_ in m:
        m_.add_particle_bonding()
        for p in m_.particles:
            if p.type.name.find('@') >= 0 and p.type.name.split('@')[0].find('H'):
                p.linker = 'head'
            elif p.type.name.find('@') >= 0 and p.type.name.split('@')[0].find('T'):
                p.linker = 'tail'
        m_.remove_linker_types()

    if s_ is None:
        s = system.replicate(m[0], 1, density=density/nmon)
    else:
        s = system.replicate(m[0], 1, s_=s_, density=density/nmon)
    print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), 1, nmon))

    for p in s.particles:
        if p.linker == 'head':
            last_head = p

        elif p.linker == 'tail':
            last_tail = p

    for m_ in m:
        if capped:
            m_.particles.remove(1)
            m_.remove_spare_bonding()
            m_.add_particle_bonding()

    s.add_particle_bonding()
    
    if traj:
        s.write_xyz('random_walk.xyz')

    temp_nmon = 1

    while True:

        m_ = m.pop(0)
        m.append(m_)
        p_ = pattern.pop(0)
        pattern.append(p_)

        if temp_nmon == 1 and p_ == 1:
            m_ = m.pop(0)
            m.append(m_)
            p_ = pattern.pop(0)
            pattern.append(p_)
        elif temp_nmon == 1:
            p_ -= 1

        for insert in range(p_):

            head = None
            tail = None

            backbone_vector = np.array([last_head.x - last_tail.x,
                                        last_head.y - last_tail.y,
                                        last_head.z - last_tail.z])

            ref_head = None
            ref_tail = None
            for p in m_.particles:
                if p.linker == 'head':
                    ref_head = p
                elif p.linker == 'tail':
                    ref_tail = p
            if ref_head and ref_tail:
                ref_backbone_vector = np.array([ref_head.x - ref_tail.x,
                                                ref_head.y - ref_tail.y,
                                                ref_head.z - ref_tail.z])
                rot_matrix = calc.find_rotation(ref_backbone_vector, backbone_vector)
                m_.rotate(around=ref_tail, rot_matrix=rot_matrix)
                translation_vector = [last_tail.x - ref_tail.x,
                                      last_tail.y - ref_tail.y,
                                      last_tail.z - ref_tail.z]
                for p in m_.particles:
                    p.x = p.x + translation_vector[0] + 3*backbone_vector[0]
                    p.y = p.y + translation_vector[1] + 3*backbone_vector[1]
                    p.z = p.z + translation_vector[2] + 3*backbone_vector[2]
            else:
                print('reference molecule has no head or tail')

            n = m_.copy()

            if capped:
                s.particles.remove(s.particles.count)
                s.remove_spare_bonding()
                s.add_particle_bonding()

            s.add(n, change_dim=False)

            s.add_particle_bonding()

            head = last_head
            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'tail':
                    tail = p

            s.make_new_bonds(head, tail, f)
            temp_nmon += 1
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), temp_nmon, nmon))

            if unwrap:
                s.unwrap()
                
            if sim is None:
                sim = lmps.Simulation(s, name='relax_%03d' % (temp_nmon), log='relax.log', **settings)
                sim.add_md(ensemble='nve', limit=limit, **settings)
                sim.add_min(**settings)
            if isinstance(sim, lmps.Simulation):
                sim.system = s
                sim.name = 'relax_%03d' % (temp_nmon)
                sim.run(np=settings.get('np'))

            if unwrap:
                s.unwrap()

            if unwrap:
                s.wrap()

            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'head':
                    last_head = p
                elif p.linker == 'tail':
                    last_tail = p

        if temp_nmon >= nmon:
            break
        
        if unwrap:
            if not s.unwrap():
                error_print('something went wrong')
                return s
    
        if traj:
            s.write_xyz('random_walk.xyz', append=True)
    
        if unwrap:
            s.wrap()
            
    for p in s.particles:
        if p not in s.molecules[p.molecule.tag].particles:
            s.molecules[p.molecule.tag].particles.add(p)

    s.write_lammps('polymer.lmps')
    s.unwrap()
    s.write_xyz('polymer.xyz')

    return s


def random_walk(m, nmon, s_=None, **kwargs):
    """pysimm.apps.random_walk.random_walk

    Builds homopolymer using random walk methodology

    Args:
        m: reference monomer :class:`~pysimm.system.System`
        nmon: total number of monomers to add to chain
        s_: :class:`~pysimm.system.System` in which to build polymer chain (None)
        extra_bonds: EXPERMINTAL, True if making ladder backbone polymer
        settings: dictionary of simulation settings
        density: density at which to build polymer (0.3)
        forcefield: :class:`~pysimm.forcefield.Forcefield` object to acquire new force field parameters
        capped: True/False if monomers are capped
        unwrap: True to unwrap final system
        traj: True to build xyz trajectory of polymer growth (True)
        limit: during MD, limit atomic displacement by this max value (LAMMPS ONLY)
        sim: :class:`~pysimm.lmps.Simulation` object for relaxation between polymer growth
    Returns:
        new polymer :class:`~pysimm.system.System`
    """
    m = m.copy()

    extra_bonds = kwargs.get('extra_bonds', False)

    settings = kwargs.get('settings', {})
    density = kwargs.get('density', 0.3)
    f = kwargs.get('forcefield')
    capped = kwargs.get('capped')
    unwrap = kwargs.get('unwrap')
    traj = kwargs.get('traj', True)
    limit = kwargs.get('limit', 0.1)
    sim = kwargs.get('sim')

    m.add_particle_bonding()

    for p in m.particles:
        if p.type.name.find('@') >= 0 and p.type.name.split('@')[0].find('H'):
            p.linker = 'head'
        elif p.type.name.find('@') >= 0 and p.type.name.split('@')[0].find('T'):
            p.linker = 'tail'

    m.remove_linker_types()

    if s_ is None:
        s = system.replicate(m, 1, density=density/nmon)
    else:
        s = system.replicate(m, 1, s_=s_, density=None)
    print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), 1, nmon))

    if traj:
        s.write_xyz('random_walk.xyz')

    if capped:
        m.particles.remove(1)
        m.remove_spare_bonding()
        m.add_particle_bonding()

    for insertion in range(nmon - 1):

        head = None
        tail = None

        backbone_vector = np.array(find_last_backbone_vector(s, m))

        for p, p_ in izip(s.particles[-1*m.particles.count:], m.particles):
                p_.x = p.x + 3*backbone_vector[0]
                p_.y = p.y + 3*backbone_vector[1]
                p_.z = p.z + 3*backbone_vector[2]

        n = m.copy()

        if capped:
            s.particles.remove(s.particles.count)
            s.remove_spare_bonding()
            s.add_particle_bonding()

        if extra_bonds:
            heads = []
            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'head':
                    heads.append(p)
        else:
            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'head':
                    head = p

        s.add(n, change_dim=False)

        s.add_particle_bonding()

        if extra_bonds:
            tails = []
            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'tail':
                    tails.append(p)
        else:
            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'tail':
                    tail = p

        for p in s.particles:
            if not p.bonded_to:
                print(p.tag)

        if head and tail:
            s.make_new_bonds(head, tail, f)
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), insertion+2, nmon))
        elif extra_bonds and len(heads) == len(tails):
            for h, t in izip(heads, tails):
                s.make_new_bonds(h, t, f)
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), insertion+2, nmon))
        else:
            print('cannot find head and tail')

        if sim is None:
            sim = lmps.Simulation(s, name='relax_%03d' % (insertion+2), log='relax.log', **settings)
            sim.add_md(ensemble='nve', limit=limit, **settings)
            sim.add_min(**settings)
        if isinstance(sim, lmps.Simulation):
            sim.system = s
            sim.name = 'relax_%03d' % (insertion+2)
            sim.run(np=settings.get('np'))

        s.unwrap()

        if traj:
            s.write_xyz('random_walk.xyz', append=True)

        if unwrap:
            s.wrap()
            
    for p in s.particles:
        if p not in s.molecules[p.molecule.tag].particles:
            s.molecules[p.molecule.tag].particles.add(p)

    s.write_lammps('polymer.lmps')
    s.unwrap()
    s.write_xyz('polymer.xyz')

    return s
