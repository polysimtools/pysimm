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


from time import strftime
from itertools import permutations, izip

import numpy as np

from pysimm import system, lmps, forcefield, calc
from pysimm import error_print


def add_bonds(s, p1, p2, f):

    if p1.molecule.particles.count < p2.molecule.particles.count:
        old_molecule_tag = p1.molecule.tag
        for p_ in p1.molecule.particles:
            p_.molecule = p2.molecule
    else:
        old_molecule_tag = p2.molecule.tag
        for p_ in p2.molecule.particles:
            p_.molecule = p1.molecule
    s.molecules.remove(old_molecule_tag)

    s.add_bond(p1, p2, f)
    for p in p1.bonded_to:
        s.add_angle(p, p1, p2, f)
        for pb in p.bonded_to:
            if pb is not p1:
                s.add_dihedral(pb, p, p1, p2, f)
    for p in p2.bonded_to:
        s.add_angle(p1, p2, p, f)
        for pb in p.bonded_to:
            if pb is not p2:
                s.add_dihedral(p1, p2, p, pb, f)
    for pb1 in p1.bonded_to:
        for pb2 in p2.bonded_to:
            s.add_dihedral(pb1, p1, p2, pb2, f)

    print('{} {}'.format(p1.tag, p2.tag))

    p1.bonded_to.add(p2)
    p2.bonded_to.add(p1)

    if s.ff_class == '2':
        for perm in permutations(p1.bonded_to, 3):
            unique = True
            for i in s.impropers:
                if i.a is not p1:
                    continue
                if set([i.b, i.c, i.d]) == set([perm[0], perm[1],
                                                perm[2]]):
                    unique = False
                    break
            if unique:
                s.add_improper(p1, perm[0], perm[1], perm[2], f)
        for perm in permutations(p2.bonded_to, 3):
            unique = True
            for i in s.impropers:
                if i.a is not p2:
                    continue
                if set([i.b, i.c, i.d]) == set([perm[0], perm[1],
                                                perm[2]]):
                    unique = False
                    break
            if unique:
                s.add_improper(p2, perm[0], perm[1], perm[2], f)


def find_last_backbone_vector(s, m):
    head_pos = [0, 0, 0]
    tail_pos = [0, 0, 0]
    for p in s.particles[-1*m.particles.count:]:
        if p.linker == 'head':
            head_pos = [p.x, p.y, p.z]
        elif p.linker == 'tail':
            tail_pos = [p.x, p.y, p.z]
    return [head_pos[0] - tail_pos[0], head_pos[1] - tail_pos[1], head_pos[2] - tail_pos[2]]


def copolymer(m, nmon, s_=None, **kwargs):

    m = [x.copy() for x in m]

    settings = kwargs.get('settings') if kwargs.get('settings') is not None else {}
    density = kwargs.get('density') or 0.3
    f = kwargs.get('forcefield')
    capped = kwargs.get('capped')
    unwrap = kwargs.get('unwrap')
    pattern = kwargs.get('pattern') or [1 for _ in range(len(m))]

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

    s.write_xyz('step_001.xyz')

    s.add_particle_bonding()

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

            add_bonds(s, head, tail, f)
            temp_nmon += 1
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), temp_nmon, nmon))

            if unwrap:
                s.unwrap()

            s.write_xyz('step_%03d.xyz' % temp_nmon)

            lmps.relax(s, dump=100, name='relax_%03d' % temp_nmon, **settings)

            lmps.minimize(s, **settings)

            if unwrap:
                s.unwrap()

            s.write_xyz('step_%03d.xyz' % temp_nmon, append=True)

            if unwrap:
                s.wrap()

            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'head':
                    last_head = p
                elif p.linker == 'tail':
                    last_tail = p

        if temp_nmon >= nmon:
            break

    s.write_lammps('copolymer.lmps')
    s.unwrap()
    s.write_xyz('copolymer.xyz')

    return s


def random_walk(m, nmon, s_=None, **kwargs):

    m = m.copy()

    extra_bonds = kwargs.get('extra_bonds') if kwargs.get('extra_bonds') is not None else False

    settings = kwargs.get('settings') if kwargs.get('settings') is not None else {}
    density = kwargs.get('density') or 0.3
    f = kwargs.get('forcefield')
    capped = kwargs.get('capped')
    unwrap = kwargs.get('unwrap')
    traj = kwargs.get('traj') if kwargs.get('traj') is not None else True

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
                p_.x = p.x + 2*backbone_vector[0]
                p_.y = p.y + 2*backbone_vector[1]
                p_.z = p.z + 2*backbone_vector[2]

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
            add_bonds(s, head, tail, f)
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), insertion+2, nmon))
        elif extra_bonds and len(heads) == len(tails):
            for h, t in izip(heads, tails):
                add_bonds(s, h, t, f)
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), insertion+2, nmon))
        else:
            print('cannot find head and tail')

        sim = lmps.Simulation(s, name='relax_%03d' % (insertion+2), log='relax.log', **settings)
        sim.add_md(ensemble='nve', limit=0.1, new_v=True, temp=300, **settings)
        sim.add_min(**settings)
        sim.run(np=settings.get('np'))

        if unwrap:
            if not s.unwrap():
                error_print('something went wrong')
                return s

        if traj:
            s.write_xyz('random_walk.xyz', append=True)

        if unwrap:
            s.wrap()

    s.write_lammps('polymer.lmps')
    s.unwrap()
    s.write_xyz('polymer.xyz')

    return s


def random_walk_multi(m, nmon, nmol=1, **kwargs):

    settings = kwargs.get('settings') if kwargs.get('settings') is not None else {}
    density = kwargs.get('density') or 0.3
    f = kwargs.get('forcefield')
    capped = kwargs.get('capped')

    distribution = kwargs.get('distribution')

    if not distribution:
        distribution = []
        for i in range(nmol):
            distribution.append(int(nmon/nmol))
        if nmon % nmol != 0:
            distribution[-1] += nmon % nmol

    m.add_particle_bonding()

    for p in m.particles:
        if p.type.name.find('@') >= 0 and p.type.name.split('@')[0].find('H'):
            p.linker = 'head'
        elif p.type.name.find('@') >= 0 and p.type.name.split('@')[0].find('T'):
            p.linker = 'tail'

    m.remove_linker_types()

    s = system.System()

    for d in distribution:

        s = system.replicate(m, 1, s_=s, density=density/nmon)
        print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), 1, nmon))

        if capped:
            m.particles.remove(1)
            m.remove_spare_bonding()
            m.add_particle_bonding()

        s.write_xyz('step_001.xyz')

        s.add_particle_bonding()

        for insertion in range(d - 1):

            head = None
            tail = None

            backbone_vector = np.array(find_last_backbone_vector(s, m))

            for p, p_ in izip(s.particles[-1*m.particles.count:], m.particles):
                p_.x = p.x + 2*backbone_vector[0]
                p_.y = p.y + 2*backbone_vector[1]
                p_.z = p.z + 2*backbone_vector[2]

            n = m.copy()

            if capped:
                s.particles.remove(s.particles.count)
                s.remove_spare_bonding()
                s.add_particle_bonding()

            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'head':
                    head = p

            s.add(n, change_dim=False)

            s.add_particle_bonding()

            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'tail':
                    tail = p

            for p in s.particles:
                if not p.bonded_to:
                    print(p.tag)

            if head and tail:
                add_bonds(s, head, tail, f)
                print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), insertion+2, nmon))
            else:
                print('cannot find head and tail')

            s.unwrap()

            s.write_xyz('step_%03d.xyz' % (insertion+2))

            lmps.relax(s, dump=100, name='relax_%03d' % (insertion+2), **settings)

            lmps.minimize(s, **settings)

            s.unwrap()

            s.write_xyz('step_%03d.xyz' % (insertion+2), append=True)

            s.wrap()

    s.write_lammps('final.lmps')
    s.unwrap()
    s.write_xyz('final.xyz')

    return s