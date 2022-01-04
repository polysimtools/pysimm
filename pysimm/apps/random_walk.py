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
import random
import sys
from time import strftime
from itertools import permutations

import numpy as np

from pysimm import system, lmps, calc
from pysimm import error_print

import math
from scipy.spatial.transform import Rotation as R
try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass


def displ_next_unit_default(m, s):
    """pysimm.apps.random_walk.displ_next_unit_default

    Default implementation of displacement of next repetitive unit in random walk method for a polymer growth

    Args:
        m: :class:`~pysimm.system.System` object -- updated
        s: :class:`~pysimm.utils.ItemContainer` of `~pysimm.system.Particle` objects
    Returns:
        empty list: general placeholder for connectivity order of linker atoms; default implementation assumes
        work with polymers with a single linkage bond, thus no connectivity order is needed
    """
    head_pos = np.zeros(3)
    tail_pos = np.zeros(3)
    hcount = 0
    tcount = 0
    bnd_lngth = 1.8

    for p in s:
        if p.linker == 'head':
            hcount += 1
            head_pos += np.array([p.x, p.y, p.z])
        elif p.linker == 'tail':
            tcount += 1
            tail_pos += np.array([p.x, p.y, p.z])

    displ = head_pos / hcount - tail_pos / tcount
    displ_dir = displ / np.linalg.norm(displ)

    for p, p_ in zip(s, m.particles):
        p_.x = p.x + displ[0] + bnd_lngth * displ_dir[0]
        p_.y = p.y + displ[1] + bnd_lngth * displ_dir[1]
        p_.z = p.z + displ[2] + bnd_lngth * displ_dir[2]
    return []


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
        geometry_rule: a pointer to a method that orients series of atoms of the next repetitive unit in random run
        settings: dictionary of simulation settings
        density: density at which to build polymer (0.3)
        forcefield: :class:`~pysimm.forcefield.Forcefield` object to acquire new force field parameters
        capped: True/False if monomers are capped
        unwrap: True to unwrap final system
        traj: True to build xyz trajectory of polymer growth (True)
        limit: during MD, limit atomic displacement by this max value (LAMMPS ONLY)
        sim: :class:`~pysimm.lmps.Simulation` object for relaxation between polymer growth
        debug: Boolean; print extra-output
    Returns:
        new polymer :class:`~pysimm.system.System`
    """
    m = m.copy()

    extra_bonds = kwargs.get('extra_bonds', False)
    displ_next_unit = kwargs.get('geometry_rule', displ_next_unit_default)

    settings = kwargs.get('settings', {})
    density = kwargs.get('density', 0.3)
    f = kwargs.get('forcefield')
    capped = kwargs.get('capped')
    unwrap = kwargs.get('unwrap')
    traj = kwargs.get('traj', True)
    limit = kwargs.get('limit', 0.1)
    sim = kwargs.get('sim')
    debug = kwargs.get('debug', False)

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

    # Remove tail-cap if it exists
    if capped:
        if __check_tags__(m, req_tags=['tail', 'tail_cap']):
            for p in m.particles:
                if p.linker == 'tail':
                    for p_ in p.bonded_to:
                        if p_.rnd_wlk_tag == 'tail_cap':
                            p.charge += p_.charge  # unite charge of tailcap into head
                            m.particles.remove(p_.tag)  # remove tailcap of monomer
                    m.remove_spare_bonding()
                    break
            m.add_particle_bonding()
        else:
            sys.exit("The capped flag is on, however, the 'tail_cap' atom is not defined")

    for insertion in range(nmon - 1):
        head = None
        tail = None

        info = displ_next_unit(m, s.particles[-1 * m.particles.count:])
        n = m.copy()

        if extra_bonds:
            heads = []
            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'head':
                    heads.append(p)
        else:
            for p in s.particles[-1*n.particles.count:]:
                if p.linker == 'head':
                    head = p

        # Remove head-cap if it exists
        if capped:
            if __check_tags__(m, req_tags=['head_cap']):
                for p_ in s.particles[-m.particles.count:]:
                    if p_.rnd_wlk_tag == 'head_cap':
                        head.charge += p_.charge  # unite charge of head_cap into tail atom
                        s.particles.remove(p_.tag)  # Removing head_cap atom from growing chain
                        s.remove_spare_bonding()
                        break
                s.add_particle_bonding()
            else:
                sys.exit("The capped flag is on, however, the 'head_cap' atom is not defined")

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

        if debug:
            for p in s.particles:
                if not p.bonded_to:
                    print(p.tag)

        if head and tail:
            s.make_new_bonds(head, tail, f)
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), insertion + 2, nmon))
        elif extra_bonds and (len(heads) == len(tails)):
            order = [(0, 0), (1, 1)]
            if len(info) == 2:
                order = [(0, info[0]), (1, info[1])]
            for elm in order:
                s.make_new_bonds(heads[elm[0]], tails[elm[1]], f)

            '''
            for h, t, ord in zip(heads, tails, extra_bonds):
                s.make_new_bonds(h, tails[ord], f)
            '''
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), insertion + 2, nmon))
            s.write_lammps('curr_progress.lmps')
        else:
            print('cannot find head and tail')

        if sim is None:
            sim = lmps.Simulation(s, name='relax_%03d' % (insertion + 2), log='relax.log', **settings)
            sim.add_md(ensemble='nve', limit=limit, **settings)
            sim.add_min(**settings)
        if isinstance(sim, lmps.Simulation):
            sim.system = s
            sim.name = 'relax_%03d' % (insertion + 2)
            sim.run(np=settings.get('np'))

        if traj:
            s.unwrap()
            s.write_xyz('random_walk.xyz', append=True)

        if unwrap:
            s.wrap()

    for p in s.particles:
        if p not in s.molecules[p.molecule.tag].particles:
            s.molecules[p.molecule.tag].particles.add(p)

    if debug:
        s.write_lammps('polymer.lmps')
        s.write_xyz('polymer.xyz')

    s.unwrap()
    return s


def find_last_tail_vector(s):
    """pysimm.apps.random_walk.find_last_tail_vector
    Finds vector defined by bond in the system between the tail atom and its capping atom.
    Requires list of particles s that formed a monomer connected on previous step of the polymerisation.

    Args:
        s: ItemContainer of :class:`~pysimm.system.Particle` objects

    Returns:
        list of vector components
    """
    result = None
    if not (__check_tags__(s, req_tags=['head', 'head_cap'])):
        print("Error: find_last_tail_vector() requires a capped monomer!")

    for p in s:
        if p.linker == 'head':
            for p_ in p.bonded_to:
                if p_.rnd_wlk_tag == 'head_cap':
                    result = [p_.x - p.x, p_.y - p.y, p_.z - p.z]
    return result


def rot_mat_about_axis(v, theta):
    """pysimm.apps.random_walk.rot_mat_about_axis
    This function returns the matrix that represents a rotation about vector v by theta degrees. Used for isotactic insertions of monomers

    Args:
        v: vector about which to rotate
        theta: degrees to rotate
    Returns:
        matrix representation of rotation
    """
    theta = theta * 2 * math.pi / 180
    r = R.from_rotvec(theta * v)
    print("Rotating vector: " + str(r.as_rotvec()))
    #as_dcm() depreciated in scipy 1.6+ #return r.as_dcm()
    return r.as_matrix()

def define_plane(a1, a2, a3):
    """pysimm.apps.random_walk.define_plane
    This function returns the mathematical constants defining a plane containing three input particles

    Args:
        a1,a2,a3: three atoms or particles
    Returns:
        np.array containing a,b,c,d that define the plane a*x + b*y + c*z + d = 0 that contains the input particles
    """
    p1 = np.array(a1.coords())
    p2 = np.array(a2.coords())
    p3 = np.array(a3.coords())
    cp = np.cross(p3 - p1, p2 - p1)
    a, b, c = cp
    d = -np.dot(cp, p3)
    return np.array([a, b, c, d])


def reflect_coords_thru_plane(atom, plane):
    """pysimm.apps.random_walk.reflect_coords_thru_plane
    This function reflects an atom through a plane, and is used for implementing syndiotactic insertions of monomers

    Args:
        atom: either an atom or an array containing x,y,z coordinates for an atom, to be reflected through the plane
        plane: np.array containing a,b,c,d that define a plane, a*x + b*y + c*z + d = 0
    Returns:
        new coordinates after reflection through plane
    """
    try:
        x1, y1, z1 = atom.coords()
    except:
        x1, y1, z1 = atom
    a, b, c, d = plane
    k = (-a * x1 - b * y1 - c * z1 - d) / float((a * a + b * b + c * c))
    x2 = a * k + x1
    y2 = b * k + y1
    z2 = c * k + z1
    x3 = 2 * x2 - x1
    y3 = 2 * y2 - y1
    z3 = 2 * z2 - z1
    # print("reflected to: " + str(atom))
    return x3, y3, z3


def scale_monomer(atom, origin, scale):
    """pysimm.apps.random_walk.scale_monomer
      This function scales the atom--origin vector. It is used by redo_monomer_insertion to scale the last monomer
      relative to its attachment point to the polymer chain

      Args:
          atom: either an atom or an array containing x,y,z coordinates for an atom, to be scaled relative to the origin
          origin: either an atom or an array containing x,y,z coordinates for where the "atom" argument should be scaled to
          scale: the factor by which the atom--origin vector should be scaled.
      Returns:
          scaled atom--origin vector
    """
    try:
        x1, y1, z1 = atom.coords()
        x0, y0, z0 = origin.coords()
    except:
        x1, y1, z1 = atom
        x0, y0, z0 = origin
    return np.array([x0 + (x1 - x0) * scale, y0 + (y1 - y0) * scale, z0 + (z1 - z0) * scale])


def redo_monomer_insertion(s, m, i):
    """pysimm.apps.random_walk.redo_monomer_insertion
    This function is called by random_walk_tacticity if the latest capped monomer insertion resulted in hardcore overlaps.
    1) The hardcore overlap is resolved by shrinking the last monomer by a factor of 0.8, iteratively, until there are no more hardcore overlaps.
    2) Then the shrunken last monomer is frozen while the rest of the polymer chain is optimized, and the last monomer is scaled in size by 1.05
    3) Cycles of contrainedOptimization and regrowth are alternated until a reasonable structure is obtained

    Args:
        s_: :class:`~pysimm.system.System` is a polymer chain in which the last monomer insertion has generated a hardcore overlap
        m: reference monomer :class:`~pysimm.system.System`. Must be a capped monomer, with headCap and tail_cap as the first and last atoms in the .mol file.
        i: number of the offending monomer, used for labelling diagnostic .xyz output files
    Returns:
        nothing; all changes to the polymer chain are written to the argument s_
    """
    for p in s.particles[-1 * m.particles.count:]:
        if p.linker == 'tail':
            tail = p
    scale_min = 0.1
    s.unwrap()
    s.set_box(padding=10)
    s.wrap()
    # shrink last monomer
    for p in s.particles[-1 * m.particles.count:]:
        p.x, p.y, p.z = scale_monomer(p, tail, scale_min)
        # now, reexpand the monomer and relax polymer, step-wise
    scale = 1
    while scale_min * scale * 1.05 < 0.91:
        print("Scaling up from %s to %s" % (str(scale_min * scale), str(scale * scale_min * 1.05)))
        scale = scale * 1.05
        for p in s.particles[-1 * m.particles.count:]:
            p.x, p.y, p.z = scale_monomer(p, tail, 1.05)
        # simulation with fixed latest monomer
        constrained_opt(s, m, "nearby")  # system-wide constrained optimization is too slow
        s.unwrap()
        s.write_xyz('bad_insertion_' + str(i) + '.xyz', append=True)
        s.wrap()
    if s.quality(tolerance=0.2) > 0:
        error_print("system is broken upon monomer reexpansion")
    # now relax the last monomer
    constrained_opt(s, m, "monomer")


def constrained_opt(s, m, active):
    """pysimm.apps.random_walk.constrained_opt
    This function is called by redo_monomer_insertion and optimizes polymer chain s while keeping the last monomer fixed.

    Args:
        s: :class:`~pysimm.system.System` is a polymer chain in which the last monomer insertion has generated a hardcore overlap
        m: reference monomer :class:`~pysimm.system.System`. Must be a capped monomer, with headCap and tail_cap as the first and last atoms in the .mol file.
    Returns:
        nothing; all changes to the polymer chain are written to the argument s_
    """
    print("Constrained Opt...")
    sim = lmps.Simulation(s, name='constrained_opt')
    total_atoms = s.particles.count
    monomer_atoms = m.particles.count
    p = s.particles[total_atoms]
    sim.add_custom("group last_monomer id " + str(total_atoms - monomer_atoms) + ":" + str(total_atoms))
    sim.add_custom(
        "group prev_two_monomers id " + str(total_atoms - 3 * monomer_atoms) + ":" + str(total_atoms - monomer_atoms))
    sim.add_custom("group non_last_monomers subtract all last_monomer")
    sim.add_custom("region insertion_area sphere {0} {1} {2} 20 side out units box".format(p.x, p.y, p.z))
    sim.add_custom("group 20_ang_away region insertion_area")
    sim.add_custom("group last_monomer_and_far union 20_ang_away last_monomer")
    if (active == "system"):
        sim.add_custom("fix freeze last_monomer setforce 0.0 0.0 0.0")
    elif (active == "monomer"):
        sim.add_custom("fix freeze non_last_monomers setforce 0.0 0.0 0.0")
    elif (active == "nearby"):
        sim.add_custom("fix freeze last_monomer_and_far setforce 0.0 0.0 0.0")
    sim.add_min(min_style="cg")
    sim.run()


def random_walk_tacticity(m, nmon, s_=None, **kwargs):
    """pysimm.apps.random_walk.random_walk_tacticity
    Builds homopolymer with controllable tacticity from capped monomer structure

    Args:
        m: reference monomer :class:`~pysimm.system.System`. Must be a capped monomer, with headCap and tail_cap
        as the first and last atoms in the .mol file.
        nmon: total number of monomers to add to chain
        s_: :class:`~pysimm.system.System` in which to build polymer chain (None)
        extra_bonds: EXPERMINTAL, True if making ladder backbone polymer
        settings: dictionary of simulation settings
        density: density at which to build polymer (0.3)
        forcefield: :class:`~pysimm.forcefield.Forcefield` object to acquire new force field parameters
        unwrap: True to unwrap final system
        debug: Boolean; print extra-output (False)
        traj: True to build xyz trajectory of polymer growth (True)
        limit: during MD, limit atomic displacement by this max value (LAMMPS ONLY)
        sim: :class:`~pysimm.lmps.Simulation` object for relaxation between polymer growth
        tacticity: float between 0 and 1.
            1 = 100% isotactic insertions
            0 = 100% syndiotactic insertions
            0.5 = equal changes of isotactic or syndiotactic insertions (i.e. atactic)
        rotation: degrees to rotate monomer per insertion
        md_spacing: how many monomer insertion steps to perform between MD relaxation steps (1)
        error_check: True/False for if monomers should be checked for hardcore overlaps after insertion
    Returns:
        new polymer :class:`~pysimm.system.System`
    """
    m = m.copy()
    extra_bonds = kwargs.get('extra_bonds', False)
    settings = kwargs.get('settings', {})
    density = kwargs.get('density', 0.3)
    f = kwargs.get('forcefield')

    unwrap = kwargs.get('unwrap')
    traj = kwargs.get('traj', True)
    debug = kwargs.get('debug', False)
    limit = kwargs.get('limit', 0.1)
    sim = kwargs.get('sim')
    tacticity = kwargs.get('tacticity', 0.5)
    if tacticity == 'atactic':
        tacticity = 0.5
    elif tacticity == 'isotactic':
        tacticity = 1
    elif tacticity == 'syndiotactic':
        tacticity = 0
    elif not (0 <= tacticity <= 1):
        sys.exit("tacticity must be a number between 0 and 1, or 'atactic' (0.5), "
                 "'isotactic' (1), or 'syndiotactic' (0)")
    tact_seq = [False] * round((nmon - 1) * tacticity) + [True] * ((nmon - 1) - round((nmon - 1) * tacticity))
    random.shuffle(tact_seq)

    rotation = kwargs.get('rotation', 0)
    md_spacing = kwargs.get('md_spacing', 1)
    error_check = kwargs.get('error_check', False)
    m.add_particle_bonding()
    if error_check:
        lmps.quick_min(m, min_style='fire')

    # Automatically redefine linkers if they have specially defined names
    for p in m.particles:
        if p.type.name.find('@') >= 0 and p.type.name.split('@')[0].find('H'):
            p.linker = 'head'
        elif p.type.name.find('@') >= 0 and p.type.name.split('@')[0].find('T'):
            p.linker = 'tail'
    m.remove_linker_types()

    # Check whether the monomer is decorated correctly
    if not __check_tags__(m.particles):
        sys.exit("random_walk:random_walk_tacticity() requires a **monomer capped with a single atom** as an input"
                 " (i.e. to model polyethylene, ethane as a monomer is required). \n"
                 "\tIn addition to 'head' and 'tail', 3 other tags should be defined: \n"
                 "\t\t(i) p.linker = 'mirror' for a particle that defines plane for iso- syndio- tactic reflection \n"
                 "\t\t(ii) p.rnd_wlk_tag = 'head_cap' and p.rnd_wlk_tag = 'tail_cap' for particles that capping head "
                 "and tail linkers correspondingly \n \t\t(see the example #13 of this distribution for details)")

    # Remove tail-cap if it exists
    for p in m.particles:
        if p.linker == 'tail':
            for p_ in p.bonded_to:
                if p_.rnd_wlk_tag == 'tail_cap':
                    p.charge += p_.charge  # unite charge of tailcap into head
                    m.particles.remove(p_.tag)  # remove tailcap of monomer
            m.remove_spare_bonding()
            break

    # Add first monomer to the output system
    if s_ is None:
        s = system.replicate(m, 1, density=density / nmon)
    else:
        s = system.replicate(m, 1, s_=s_, density=None)
    print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), 1, nmon))
    if traj:
        s.write_xyz('random_walk.xyz')
    s.add_particle_bonding()

    # Main polymerisation loop
    for insertion in range(nmon - 1):
        n = m.copy()
        head = None
        tail = None
        mirror_atom = None
        for p in n.particles:
            if p.linker == 'head':
                head = p
            elif p.linker == 'tail':
                tail = p
            elif p.linker == 'mirror':
                mirror_atom = p
        backbone_vector = np.array(find_last_backbone_vector(s, m))
        tail_vector = np.array(find_last_tail_vector(s.particles[-n.particles.count:]))

        for p, p_ in zip(s.particles[-1 * n.particles.count:], n.particles):  # translate monomer
            a = 1.1  # coefficient of displacement of a new monomer along the head--tail direction
            b = 2.4  # coefficient of displacement of a new monomer along the head--headcap direction
            p_.x = p.x + a * backbone_vector[0] + b * tail_vector[0]
            p_.y = p.y + a * backbone_vector[1] + b * tail_vector[1]
            p_.z = p.z + a * backbone_vector[2] + b * tail_vector[2]

        if tact_seq[insertion]:  # if syndiotactic insertion, reflect monomer
            print("syndiotactic insertion...")
            mirrorPlane = define_plane(head, tail, mirror_atom)
            for p in n.particles:
                p.x, p.y, p.z = reflect_coords_thru_plane([p.x, p.y, p.z], mirrorPlane)

        else:  # else isotatic insertion, rotate monomer if necessary
            print("isotatic insertion...")
            if rotation != 0:  # rotate monomer, if necessary
                rot_mat = rot_mat_about_axis(backbone_vector, rotation)
                n.rotate(around=head, rot_matrix=rot_mat)

        for p_ in s.particles[-n.particles.count:]:
            if p_.rnd_wlk_tag == 'head_cap':
                head.charge += p_.charge  # unite charge of head_cap into tail atom
                s.particles.remove(p_.tag)  # Removing head_cap atom from growing chain
                s.remove_spare_bonding()
                break

        if extra_bonds:
            heads = []
            for p in s.particles[-n.particles.count:]:
                if p.linker == 'head':
                    heads.append(p)
        else:
            for p in s.particles[-n.particles.count:]:
                if p.linker == 'head':
                    head = p

        s.add(n, change_dim=False)
        s.add_particle_bonding()
        if extra_bonds:
            tails = []
            for p in s.particles[-n.particles.count:]:
                if p.linker == 'tail':
                    tails.append(p)
        else:
            for p in s.particles[-n.particles.count:]:
                if p.linker == 'tail':
                    tail = p

        if debug:
            for p in s.particles:
                if not p.bonded_to:
                    print(p.tag)

        if head and tail:
            s.make_new_bonds(head, tail, f)
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), insertion + 2, nmon))
        elif extra_bonds and len(heads) == len(tails):
            for h, t in zip(heads, tails):
                s.make_new_bonds(h, t, f)
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), insertion + 2, nmon))
        else:
            print('cannot find head and tail')
        if sim is None:
            sim = lmps.Simulation(s, name='relax_%03d' % (insertion + 2), log='relax.log', **settings)
            if (insertion + 2) % md_spacing == 0:
                sim.add_md(ensemble='nve', limit=limit, **settings)
            # sim.add_min(**settings)

        if isinstance(sim, lmps.Simulation):
            s_ = s.copy()
            sim.system = s
            sim.name = 'relax_%03d' % (insertion + 2)
            sim.run(np=settings.get('np'))
            energy = lmps.energy(s)
            print("LAMMPS Energy = " + str(energy))
            print("LAMMPS Energy/#ofAtoms = " + str(energy / s.particles.count))
            if error_check == True:  # check for hardcore overlap
                print("checking for hardcore overlap")
                if s.quality(tolerance=0.3) > 0:
                    print("Found bad quality monomer insertion. Redoing last insertion...")
                    s.unwrap()
                    s.write_xyz('bad_insertion_' + str(insertion + 2) + '.xyz')
                    s.wrap()
                    redo_monomer_insertion(s_, n, insertion + 2)
                    s = s_.copy()
        if traj:
            s.unwrap()
            s.write_xyz('random_walk.xyz', append=True)

    # Removing the very last 'head_cap' at the end of the chain
    for p_ in s.particles[-n.particles.count:]:
        if p_.rnd_wlk_tag == 'head_cap':
            head.charge += p_.charge  # unite charge of head_cap into tail atom
            s.particles.remove(p_.tag)  # Removing head_cap atom from growing chain
            s.remove_spare_bonding()
    # Syncronizing molecule representation with particles ItemContainer representation for the chain
    s.objectify()

    if debug:
        s.write_lammps('polymer.lmps')
        s.write_xyz('polymer.xyz')

    s.unwrap()
    return s


def __check_tags__(m, **kwargs):
    """pysimm.apps.random_walk.__check_tags__
        private method to assert the polymerisation-related decorators assigned to the system 'm' that represents the
        next repetitive unit
    """
    tags = [p.linker for p in m] + [p.rnd_wlk_tag for p in m]
    req_tags = kwargs.get('req_tags', ['head', 'tail', 'mirror', 'head_cap', 'tail_cap'])
    tmp = True
    for tg in req_tags:
        tmp *= (tg in tags)
    return bool(tmp)


def check_tacticity(s, char_idxs, mon_len):
    """pysimm.apps.random_walk.check_tacticity
        Method evaluates the local geometry of the polymer :class:`~pysimm.system.System`.
        correct input includes
        Args:
            char_idxs (list of int): characteristic indexes that define the structure of repetetive unit of the monomer.
            It is supposed to have 4 elements which define index of (1) first atom in backbone; (2) second atom in the
            backbone; (3) closest to backbone atom on the fisrt side chain; (4) closest to backbone atom on the second
            side chain
            mon_len (int): number of atoms in uncapped rep. unit

        Note: currentely it is assumed that polymerisation does not change local indexing so indexes of corresponding
        characteristic atoms of the chain can be found by adding a number multiple of mon_len

        Returns:
            angles (list of float): angles (in deg) between corresponding pairs of backbone vector (1-2) and normal to
            the plane produced by to side chains (2-3 x 2-4). Those vectors can be either on one half-space of (3-2-4)
            plane, so the angle will be >90 (deg) or on the opposite half-spaces of the plane, so the angle <90 (deg).
            orientations (list of boolean): sequence that tracks local geometry of a chain: records True if two
            consecutive rep.units form a meso dyad, and False if they form a racemo dyad
    """

    offset = 0
    # backbone vector pnt-1
    tails = [s.particles[offset + char_idxs[0]]] + \
            [s.particles[i] for i in range(offset + char_idxs[0] + mon_len, len(s.particles), mon_len)]
    # backbone vector pnt-2
    heads = [s.particles[offset + char_idxs[1]]] + \
            [s.particles[i] for i in range(offset + char_idxs[1] + mon_len, len(s.particles), mon_len)]
    # start side chain-1
    sides = [s.particles[offset + char_idxs[2]]] + \
            [s.particles[i] for i in range(offset + char_idxs[2] + mon_len, len(s.particles), mon_len)]
    # start side chain-2
    methyls = [s.particles[offset + char_idxs[3]]] + \
              [s.particles[i] for i in range(offset + char_idxs[3] + mon_len, len(s.particles), mon_len)]

    angles = []
    for h, t, s, m in zip(heads, tails, sides, methyls):
        HT = np.array([h.x - t.x, h.y - t.y, h.z - t.z])
        HT /= np.linalg.norm(HT)

        Hmethyl = np.array([h.x - m.x, h.y - m.y, h.z - m.z])
        Hmethyl /= np.linalg.norm(Hmethyl)

        Hside = np.array([h.x - s.x, h.y - s.y, h.z - s.z])
        Hside /= np.linalg.norm(Hside)

        side_X_methyl = np.cross(Hside, Hmethyl)

        side_dot_methyl = np.dot(Hside, Hmethyl)
        side_theta_methyl = np.arccos(side_dot_methyl)

        side_X_methyl /= np.sin(side_theta_methyl)

        cos_theta = np.dot(side_X_methyl, HT)
        angles.append(np.arccos(cos_theta) / np.pi * 180)

        tmp = np.array([2 * int(a >= 90) - 1 for a in angles])

    return angles, [t > 0 for t in tmp[:-1] * tmp[1:]]

