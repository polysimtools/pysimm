# ******************************************************************************
# pysimm.apps.polymerize module
# ******************************************************************************
#
# Polymatic algorithm written using pysimm tools
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

import os
import sys
from time import strftime
from itertools import permutations

from pysimm import system, lmps
from pysimm.utils import ItemContainer

rappture = True
try:
    import Rappture
except ImportError:
    rappture = False


head_linkers = ItemContainer()
tail_linkers = ItemContainer()
settings = None


def initialize(s, polym_settings):
    global head_linkers, tail_linkers, settings

    settings = polym_settings

    print('finding linker atoms in system and adding artifical charges...')

    s.add_particle_bonding()

    for p in s.particles:
        if p.type.name == settings.polym.lo1:
            p.linker = True
            head_linkers.add(p)
            if settings.polym.charge:
                p.charge += settings.polym.charge
        elif p.type.name == settings.polym.lo2:
            p.linker = True
            tail_linkers.add(p)
            if settings.polym.charge:
                    p.charge -= settings.polym.charge

    s.remove_linker_types()


def finalize(s):
    global head_linkers, tail_linkers, settings

    print('renaming unbonded linker atoms and resetting charges')

    for p in head_linkers:
        old_head = s.particle_types.get(settings.polym.lo1)
        if old_head:
            p.type = old_head[0]
        if settings.polym.charge:
            p.charge -= settings.polym.charge
    for p in tail_linkers:
        old_tail = s.particle_types.get(settings.polym.lo2)
        if old_tail:
            p.type = old_tail[0]
        if settings.polym.charge:
            p.charge += settings.polym.charge

    s.unwrap()
    s.write_lammps('final.lmps')


def check_for_bonds(s):
    global head_linkers, tail_linkers, settings

    for p1 in head_linkers:
        for p2 in tail_linkers:
            if p1.molecule is p2.molecule:
                continue
            if (settings.polym.cutoff and
                    s.distance(p1, p2) < settings.polym.cutoff):
                if settings.polym.extra_bond:
                    other_head = other_tail = None
                    for p in p1.molecule.particles:
                        if p in head_linkers and p is not p1:
                            other_head = p
                            break
                    for p in p2.molecule.particles:
                        if p in tail_linkers and p is not p2:
                            other_tail = p
                            break

                    if (other_head and other_tail and
                            (s.distance(other_head, other_tail) <
                             settings.polym.cutoff)):

                        add_bonds(s, p1, p2)

                        head_linkers.remove(p1.tag, update=False)
                        if settings.polym.charge:
                            p1.charge -= settings.polym.charge

                        tail_linkers.remove(p2.tag, update=False)
                        if settings.polym.charge:
                            p2.charge += settings.polym.charge

                        add_bonds(s, other_head, other_tail)

                        head_linkers.remove(other_head.tag, update=False)
                        if settings.polym.charge:
                            other_head.charge -= settings.polym.charge
                        tail_linkers.remove(other_tail.tag, update=False)
                        if settings.polym.charge:
                            other_tail.charge += settings.polym.charge
                        return True
                elif settings.polym.max_chain_length:
                    if (p1.molecule.nmon and p2.molecule.nmon and
                            (p1.molecule.nmon + p2.molecule.nmon) <=
                            settings.polym.max_chain_length):
                        add_bonds(s, p1, p2)
                        head_linkers.remove(p1.tag, update=False)
                        tail_linkers.remove(p2.tag, update=False)
                        if settings.polym.charge:
                            p1.charge -= settings.polym.charge
                            p2.charge += settings.polym.charge
                        return True
                    elif (not p1.molecule.nmon and p2.molecule.nmon and
                            (p2.molecule.nmon + 1) <=
                            settings.polym.max_chain_length):
                        add_bonds(s, p1, p2)
                        head_linkers.remove(p1.tag, update=False)
                        tail_linkers.remove(p2.tag, update=False)
                        if settings.polym.charge:
                            p1.charge -= settings.polym.charge
                            p2.charge += settings.polym.charge
                        return True
                    elif (p1.molecule.nmon and not p2.molecule.nmon and
                            (p1.molecule.nmon + 1) <=
                            settings.polym.max_chain_length):
                        add_bonds(s, p1, p2)
                        head_linkers.remove(p1.tag, update=False)
                        tail_linkers.remove(p2.tag, update=False)
                        if settings.polym.charge:
                            p1.charge -= settings.polym.charge
                            p2.charge += settings.polym.charge
                        return True
                    elif not p1.molecule.nmon and not p2.molecule.nmon:
                        add_bonds(s, p1, p2)
                        head_linkers.remove(p1.tag, update=False)
                        tail_linkers.remove(p2.tag, update=False)
                        if settings.polym.charge:
                            p1.charge -= settings.polym.charge
                            p2.charge += settings.polym.charge
                        return True
                else:
                    add_bonds(s, p1, p2)
                    head_linkers.remove(p1.tag, update=False)
                    tail_linkers.remove(p2.tag, update=False)
                    if settings.polym.charge:
                        p1.charge -= settings.polym.charge
                        p2.charge += settings.polym.charge
                    return True

    return False


def add_bonds(s, p1, p2):
    global settings

    s.check_items()

    f = settings.polym.forcefield

    m1 = p1.molecule
    m2 = p2.molecule

    if m1.particles.count > m2.particles.count:
        large_m, small_m = m1, m2
    else:
        large_m, small_m = m2, m1

    if m1 is not m2:
        large_m.polymer = True
        if small_m.polymer and small_m.nmon and large_m.nmon:
            large_m.nmon += small_m.nmon
        elif large_m.nmon:
            large_m.nmon += 1
        else:
            large_m.nmon = 2

    s.add_bond(p1, p2, f)
    if m1 is not m2:
        large_m.bonds.add(s.bonds[-1])
    for p in p1.bonded_to:
        s.add_angle(p, p1, p2, f)
        if m1 is not m2:
            large_m.angles.add(s.angles[-1])
        for pb in p.bonded_to:
            if pb is not p1:
                s.add_dihedral(pb, p, p1, p2, f)
                if m1 is not m2:
                    large_m.dihedrals.add(s.dihedrals[-1])
    for p in p2.bonded_to:
        s.add_angle(p1, p2, p, f)
        if m1 is not m2:
            large_m.angles.add(s.angles[-1])
        for pb in p.bonded_to:
            if pb is not p2:
                s.add_dihedral(p1, p2, p, pb, f)
                if m1 is not m2:
                    large_m.dihedrals.add(s.dihedrals[-1])
    for pb1 in p1.bonded_to:
        for pb2 in p2.bonded_to:
            s.add_dihedral(pb1, p1, p2, pb2, f)
            if m1 is not m2:
                large_m.dihedrals.add(s.dihedrals[-1])

    p1.bonded_to.append(p2)
    p2.bonded_to.append(p1)

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

    if m1 is not m2:
        for p in small_m.particles:
            p.molecule = large_m
            large_m.particles.add(p)
        for b in small_m.bonds:
            large_m.bonds.add(b)
        for a in small_m.angles:
            large_m.angles.add(a)
        for d in small_m.dihedrals:
            large_m.dihedrals.add(d)
        for i in small_m.impropers:
            large_m.impropers.add(i)

    if m1 is not m2:
        s.molecules.remove(small_m.tag)


def polymerize(s, polym_settings):
    global head_linkers, tail_linkers, settings

    if rappture:
        Rappture.Utils.progress(0, 'Initializing Polymatic...')
    bonds = 0

    initialize(s, polym_settings)

    if rappture:
        Rappture.Utils.progress(0, '%s/%s bonds made: '
                                   'optimizing initial structure...'
                                % (bonds, settings.polym.max_bonds))

    os.mkdir('logs')

    lmps_min(s, 'initial optimization')

    while bonds < settings.polym.max_bonds:
        attempt = 0
        result = check_for_bonds(s)
        while not result and attempt < settings.polym.max_md:
            attempt += 1
            if rappture:
                Rappture.Utils.progress(int(float(bonds)/settings.
                                            polym.max_bonds *
                                            100),
                                        '%s/%s bonds made: attempt #%s to make '
                                        'new bond'
                                        % (bonds, settings.polym.max_bonds,
                                           attempt))
            lmps_step_md(s, bonds, attempt)
            result = check_for_bonds(s)
        if not result:
            break

        bonds += 1
        if rappture:
            Rappture.Utils.progress(int(float(bonds)/settings.polym.max_bonds
                                        * 100),
                                    '%s/%s bonds made: '
                                    'optimizing newly formed bond'
                                    % (bonds, settings.polym.max_bonds))
        print('%s: bond %s made successfully' % (strftime('%H:%M:%S'), bonds))
        sys.stdout.flush()
        monomers = 0
        polymers = []
        for m in s.molecules:
            if not m.nmon:
                monomers += 1
            else:
                polymers.append(m.nmon)
        print('%s: polymer distribution:' % strftime('%H:%M:%S'))
        for n in reversed(sorted(list(set(polymers)))):
            print('\t%s chains with %s monomers' % (polymers.count(n), n))
        print('\t%s free monomers' % monomers)
        lmps_min(s, 'bond %s optimization' % bonds)
        s.write_lammps('step_%03d.lmps' % bonds)

        if (bonds % settings.polym.cycle == 0 and
                (bonds / settings.polym.cycle) % settings.polym.npt_freq == 0):
            if rappture:
                Rappture.Utils.progress(int(float(bonds)/settings.
                                            polym.max_bonds
                                            * 100),
                                        '%s/%s bonds made: '
                                        'performing npt cycle md'
                                        % (bonds, settings.polym.max_bonds))
            lmps_cycle_npt_md(s, bonds)

        elif bonds % settings.polym.cycle == 0:
            if rappture:
                Rappture.Utils.progress(int(float(bonds)/settings.
                                            polym.max_bonds
                                            * 100),
                                        '%s/%s bonds made: '
                                        'performing nvt cycle md'
                                        % (bonds, settings.polym.max_bonds))
            lmps_cycle_nvt_md(s, bonds)

    if rappture:
        Rappture.Utils.progress(99, 'Finalizing Polymatic')

    finalize(s)


def lmps_min(s, name):
    global settings

    if settings.polym.min.cluster:
        nanohub = {'cores': int(settings.polym.min.nanohub_cores),
                   'walltime': int(settings.polym.min.nanohub_walltime)}
    else:
        nanohub = {}

    if settings.polym.min.user_input:
        lmps.run(s, lammps_in=settings.polym.min.min_in,
                 name='initial optimization',
                 print_to_screen=False, nanohub=nanohub)
    else:
        lmps.minimize(s, name=name,
                      cutoff=settings.polym.min.nb_cutoff,
                      sd_etol=settings.polym.min.sd_etol,
                      sd_ftol=settings.polym.min.sd_ftol,
                      sd_maxiter=settings.polym.min.sd_maxiter,
                      sd_maxeval=settings.polym.min.sd_maxeval,
                      cg_etol=settings.polym.min.cg_etol,
                      cg_ftol=settings.polym.min.cg_ftol,
                      cg_maxiter=settings.polym.min.cg_maxiter,
                      cg_maxeval=settings.polym.min.cg_maxeval,
                      log='logs/%s' % '_'.join(name.split()),
                      np=settings.np,
                      nanohub=nanohub)


def lmps_step_md(s, bonds, attempt):
    global settings

    if settings.polym.step.cluster:
        nanohub = {'cores': int(settings.polym.step.nanohub_cores),
                   'walltime': int(settings.polym.step.nanohub_walltime)}
    else:
        nanohub = {}

    if settings.polym.step.user_input:
        lmps.run(s, lammps_in=settings.polym.step.step_in,
                 name='bond %s attempt #%d' % (bonds + 1, attempt),
                 print_to_screen=False, nanohub=nanohub)
    else:
        lmps.md(s, name='bond %s: attempt #%d' % (bonds + 1, attempt),
                ensemble='nvt',
                cutoff=settings.polym.step.nb_cutoff,
                temp=settings.polym.step.temp,
                new_v=True,
                length=settings.polym.step.length,
                log='logs/step_%03d_%03d' % (bonds, attempt),
                dump=settings.polym.dump,
                dump_name=settings.polym.dump_name,
                dump_append=settings.polym.dump_append,
                np=settings.np,
                nanohub=nanohub)


def lmps_cycle_nvt_md(s, bonds):
    global settings

    if settings.polym.cycle_nvt.cluster:
        nanohub = {'cores': int(settings.polym.cycle_nvt.nanohub_cores),
                   'walltime': int(settings.polym.cycle_nvt.nanohub_walltime)}
    else:
        nanohub = {}

    if settings.polym.cycle_nvt.user_input:
        lmps.run(s, lammps_in=settings.polym.cycle_nvt.step_in,
                 name='bond %d cycle nvt' % bonds,
                 print_to_screen=False, nanohub=nanohub)
    else:
        lmps.md(s, name='bond %d cycle nvt' % bonds,
                ensemble='nvt',
                cutoff=settings.polym.cycle_nvt.nb_cutoff,
                temp=settings.polym.cycle_nvt.temp,
                new_v=True,
                length=settings.polym.cycle_nvt.length,
                log='logs/cycle_nvt_%03d' % bonds,
                dump=settings.polym.dump,
                dump_name=settings.polym.dump_name,
                dump_append=settings.polym.dump_append,
                np=settings.np,
                nanohub=nanohub)


def lmps_cycle_npt_md(s, bonds):
    global settings

    if settings.polym.cycle_npt.cluster:
        nanohub = {'cores': int(settings.polym.cycle_npt.nanohub_cores),
                   'walltime': int(settings.polym.cycle_npt.nanohub_walltime)}
    else:
        nanohub = {}

    if settings.polym.cycle_npt.user_input:
        lmps.run(s, lammps_in=settings.polym.cycle_npt.step_in,
                 name='bond %d cycle npt' % bonds,
                 print_to_screen=False, nanohub=nanohub)
    else:
        lmps.md(s, name='bond %d cycle npt' % bonds,
                ensemble='npt',
                cutoff=settings.polym.cycle_npt.nb_cutoff,
                temp=settings.polym.cycle_npt.temp,
                new_v=True,
                pressure=settings.polym.cycle_npt.pressure,
                length=settings.polym.cycle_npt.length,
                log='logs/cycle_npt_%03d' % bonds,
                dump=settings.polym.dump,
                dump_name=settings.polym.dump_name,
                dump_append=settings.polym.dump_append,
                np=settings.np,
                nanohub=nanohub)


def lmps_final_md(s):
    global settings

    if settings.polym.step.cluster:
        nanohub = {'cores': int(settings.polym.step.nanohub_cores),
                   'walltime': int(settings.polym.step.nanohub_walltime)}
    else:
        nanohub = {}

    lmps.md(s, name='final',
            ensemble='nvt',
            cutoff=settings.polym.step.nb_cutoff,
            temp=settings.polym.step.temp,
            new_v=True,
            length=settings.polym.step.length,
            log='logs/final',
            dump=False,
            np=settings.np,
            nanohub=nanohub)
