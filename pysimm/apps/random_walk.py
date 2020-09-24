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
# Copyright (c) 2016 Michael E. Fortunato, Coray M. Colina, Sibo Lin
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
from itertools import permutations

import numpy as np

from pysimm import system, lmps, forcefield, calc
from pysimm import error_print

import math
from scipy.spatial.transform import Rotation as R
try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass

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

        for p, p_ in zip(s.particles[-1*m.particles.count:], m.particles):
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
            for h, t in zip(heads, tails):
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

def find_last_tail_vector(s, m, **kwargs):
    """pysimm.apps.random_walk.find_last_tail_vector
    Finds vector defined by bond in the system between the tail atom and its capping atom. Requires current system s and reference monomer m. The monomer must be capped.
    
    Args:
        s: :class:`~pysimm.system.System` object
        m: :class:`~pysimm.system.System` object
    Returns:
        list of vector components
    """
    capped = kwargs.get('isCapped',False)
    if not capped:
        print("Error: find_last_tail_vector() requires a capped monomer!")
        return
    headCap = s.particles[-1]
    for p in s.particles[-1*m.particles.count:]:
        if p.linker == 'head':
            #print("found headAtom @ " + str(p.coords()))
            #print("found headCap  @ " + str(headCap.coords())) 
            return [headCap.x - p.x, headCap.y - p.y, headCap.z - p.z]
    print("no head atom found in system! check your monomer definition!")
    return

def rotMatAboutAxis(v,theta):
    """pysimm.apps.random_walk.rotMatAboutAxis
    This function returns the matrix that represents a rotation about vector v by theta degrees. Used for isotactic insertions of monomers
    
    Args:
        v: vector about which to rotate
        theta: degrees to rotate
    Returns:
        matrix representation of rotation
    """
    theta = theta * 2* math.pi / 180
    r = R.from_rotvec(theta * v)
    print("Rotating vector: " + str(r.as_rotvec()))
    return r.as_dcm()

def definePlane(a1,a2,a3):
    """pysimm.apps.random_walk.definePlane
    This function returns the mathematical constants defining a plane containing three input particles
    
    Args:
        a1,a2,a3: three atoms or particles
    Returns:
        np.array containing a,b,c,d that define the plane a*x + b*y + c*z + d = 0 that contains the input particles
    """
    p1=np.array(a1.coords())
    p2=np.array(a2.coords())
    p3=np.array(a3.coords())
    v1 = p3-p1
    v2 = p2-p1
    cp = np.cross(v1,v2)
    a,b,c = cp
    d = -np.dot(cp,p3)
    #print str(cp) + str(d)
    return np.array([a,b,c,d])
  
def reflectCoordsThruPlane(atom,plane):
    """pysimm.apps.random_walk.reflectCoordsThruPlane
    This function reflects an atom through a plane, and is used for implementing syndiotactic insertions of monomers
    
    Args:
        atom: either an atom or an array containing x,y,z coordinates for an atom, to be reflected through the plane
        plane: np.array containing a,b,c,d that define a plane, a*x + b*y + c*z + d = 0
    Returns:
        new coordinates after reflection through plane
    """
    try:
        x1,y1,z1 = atom.coords()
    except:
        x1,y1,z1 = atom
    a,b,c,d = plane
    k = (-a*x1 - b*y1 -c*z1 - d)/float((a*a + b*b + c*c))
    x2 = a*k + x1
    y2 = b*k + y1
    z2 = c*k + z1
    x3 = 2*x2 - x1
    y3 = 2*y2 - y1
    z3 = 2*z2 - z1
    #print("reflected to: " + str(atom))
    return x3,y3,z3
    
def scaleMonomer(atom,origin,scale):
  """pysimm.apps.random_walk.scaleMonomer
    This function scales the atom--origin vector. It is used by redoMonomerInsertion to scale the last monomer relative to its attachment point to the polymer chain
    
    Args:
        atom: either an atom or an array containing x,y,z coordinates for an atom, to be scaled relative to the origin
        origin: either an atom or an array containing x,y,z coordinates for where the "atom" argument should be scaled to
        scale: the factor by which the atom--origin vector should be scaled. 
    Returns:
        scaled atom--origin vector
  """
  try:
    x1,y1,z1 = atom.coords()
    x0,y0,z0 = origin.coords()
  except:
    x1,y1,z1 = atom
    x0,y0,z0 = origin
  return np.array([x0+(x1-x0)*scale,y0+(y1-y0)*scale,z0+(z1-z0)*scale])

def redoMonomerInsertion(s_,m,i):
    """pysimm.apps.random_walk.redoMonomerInsertion
    This function is called by random_walk_tacticity if the latest capped monomer insertion resulted in hardcore overlaps. 
    1) The hardcore overlap is resolved by shrinking the last monomer by a factor of 0.8, iteratively, until there are no more hardcore overlaps.
    2) Then the shrunken last monomer is frozen while the rest of the polymer chain is optimized, and the last monomer is scaled in size by 1.05
    3) Cycles of contrainedOptimization and regrowth are alternated until a reasonable structure is obtained
    
    Args:
        s_: :class:`~pysimm.system.System` is a polymer chain in which the last monomer insertion has generated a hardcore overlap
        m: reference monomer :class:`~pysimm.system.System`. Must be a capped monomer, with headCap and tailCap as the first and last atoms in the .mol file.
        i: number of the offending monomer, used for labelling diagnostic .xyz output files
    Returns:
        nothing; all changes to the polymer chain are written to the argument s_ 
    """
    s = s_.copy()
    for p in s.particles[-1*m.particles.count:]:
        if p.linker == 'tail':
            tail = p
    scale = 1
    successfulInsertion = False
    s.unwrap()
    s.set_box(padding=10)
    s.wrap()
    while successfulInsertion == False:
        scale = scale * 0.8
        print("scale = " + str(scale))
        for p in s.particles[-1*n.particles.count:]:
            p.x,p.y,p.z = scaleMonomer(p,tail,scale)
        #simulation with fixed latest monomer
        constrainedOpt(s,m)
        s.unwrap()
        s.write_xyz('bad_insertion_' + str(i) + '.xyz', append=True)
        s.wrap()
        quality = s.quality(tolerance=0.3)
        print("Quality: " + str(quality))
        if quality == 0:
            successfulInsertion = True
    #now, reexpand the monomer
    scaleMin = scale
    scale = 1
    increments=10
    s_ = s.copy()
    while scaleMin*scale*1.05 < 0.9:
        scale = scale*1.05
        print("Scaling up from %s to %s" % (str(scaleMin), str(scale*scaleMin)))
        for p in s.particles[-1*n.particles.count:]:
            p.x,p.y,p.z = scaleMonomer(p,tail,scale)
        #simulation with fixed latest monomer
        constrainedOpt(s,m)
        s.unwrap()
        s.write_xyz('bad_insertion_' + str(i) + '.xyz', append=True)
        s.wrap()
        if s.quality(tolerance=0.3) > 0:
            print("system is broken upon monomer reexpansion")
            break
    s_ = s.copy()
    
def constrainedOpt(s,m):
    """pysimm.apps.random_walk.constrainedOpt
    This function is called by redoMonomerInsertion and optimizes polymer chain s while keeping the last monomer fixed.
    
    Args:
        s: :class:`~pysimm.system.System` is a polymer chain in which the last monomer insertion has generated a hardcore overlap
        m: reference monomer :class:`~pysimm.system.System`. Must be a capped monomer, with headCap and tailCap as the first and last atoms in the .mol file.
    Returns:
        nothing; all changes to the polymer chain are written to the argument s_ 
    """
    print("Constrained Opt...")
    sim=lmps.Simulation(s,name='constrainedOpt')
    totalAtoms = s.particles.count
    monomerAtoms = m.particles.count
    sim.add_custom("group lastMonomer id " + str(totalAtoms-monomerAtoms) + ":" + str(totalAtoms))
    sim.add_custom("group otherMonomers subtract all lastMonomer")
    sim.add_custom("fix freeze lastMonomer setforce 0.0 0.0 0.0")
    sim.add_min()
    sim.run()
    
def random_walk_tacticity(m, nmon, s_=None, **kwargs):
    """pysimm.apps.random_walk.random_walk_tacticity
    Builds homopolymer with controllable tacticity from capped monomer structure
    
    Args:
        m: reference monomer :class:`~pysimm.system.System`. Must be a capped monomer, with headCap and tailCap as the first and last atoms in the .mol file.
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
        tacticity: float between 0 and 1. 
            1 = 100% isotactic insertions
            0 = 100% syndiotactic insertions
            0.5 = equal changes of isotactic or syndiotactic insertions (i.e. atactic)
        rotation: degrees to rotate monomer per insertion
        md_spacing: how many monomer insertion steps to perform between MD relaxation steps (1)
        errorCheck: True/False for if monomers should be checked for hardcore overlaps after insertion
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
    tacticity = kwargs.get('tacticity',0.5)
    if tacticity == 'atactic':
        tacticity = 0.5
    elif tacticity == 'isotactic':
        tacticity = 1
    elif tacticity == 'syndiotactic':
        tacticity = 0
    elif not ( 0 <= tacticity <= 1):
        print("tacticity must be a number between 0 and 1, or 'atactic' (0.5), 'isotactic' (1), or 'syndiotactic' (0)")
    rotation = kwargs.get('rotation',0)
    md_spacing = kwargs.get('md_spacing',1)
    errorCheck = kwargs.get('errorCheck',False)
    m.add_particle_bonding()
    if errorCheck:
        lmps.quick_min(m, min_style='fire')
        E_monomer = lmps.energy(m)
        E_system = E_monomer
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
    if not capped:
        print("random_walk_tacticity() requires capped monomers (i.e. to model polyethylene, the monomer should be ethane, where the first and last atoms in the .mol file are hydrogen atoms)")
        return
    
    for p in m.particles:
        if p.linker == 'head':
            p.charge += m.particles[1].charge #unite charge of headcap into head
            break
    m.particles.remove(1) #remove headcap of monomer
    m.remove_spare_bonding()
    m.add_particle_bonding()
    
    for insertion in range(nmon - 1):
        n = m.copy()
        head = None
        tail = None
        mirrorAtom = None
        tailCap = n.particles[-1]
        for p in n.particles:
            if p.linker == 'head':
                head = p
            elif p.linker == 'tail':
                tail = p
            elif p.linker == 'mirror':
                mirrorAtom = p
        backbone_vector = np.array(find_last_backbone_vector(s, m))
        tail_vector = np.array(find_last_tail_vector(s,m,isCapped=capped))
                
        for p, p_ in zip(s.particles[-1*n.particles.count:], n.particles): #translate monomer
            a=1
            b=1.4 #scale to convert C-H bond to C-C bond
            p_.x = p.x + a*backbone_vector[0] + b*tail_vector[0]
            p_.y = p.y + a*backbone_vector[1] + b*tail_vector[1]
            p_.z = p.z + a*backbone_vector[2] + b*tail_vector[2]
        if np.random.rand() > tacticity: #if syndiotactic insertion, reflect monomer
            print("syndiotactic insertion...")
            mirrorPlane = definePlane(head, tail, mirrorAtom)
            for p in n.particles:
                p.x,p.y,p.z = reflectCoordsThruPlane([p.x,p.y,p.z],mirrorPlane)
            
        else: #else isotatic insertion, rotate monomer if necessary
            print("isotatic insertion...")
            if rotation != 0: #rotate monomer, if necessary                                               
                rotMat = rotMatAboutAxis(backbone_vector,rotation)
                n.rotate(around=head,rot_matrix=rotMat)
        
        tail.charge += s.particles[-1].charge #unite charge of tailcap into tail atom
        s.particles.remove(s.particles.count) #Removing tailCap atom from growing chain
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
        
        sCopy = s.copy()
        
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
            for h, t in zip(heads, tails):
                s.make_new_bonds(h, t, f)
            print('%s: %s/%s monomers added' % (strftime('%H:%M:%S'), insertion+2, nmon))
        else:
            print('cannot find head and tail')
        if sim is None:
            sim = lmps.Simulation(s, name='relax_%03d' % (insertion+2), log='relax.log', **settings)
            if (insertion+2)%md_spacing == 0:
                sim.add_md(ensemble='nve', limit=limit, **settings)
            sim.add_min(**settings)
        if isinstance(sim, lmps.Simulation):
            sim.system = s
            sim.name = 'relax_%03d' % (insertion+2)
            sim.run(np=settings.get('np'))
            energy = lmps.energy(s)
            print("LAMMPS Energy = " + str(energy))
            print("LAMMPS Energy/#ofAtoms = " + str(energy/s.particles.count))
            if errorCheck == True: #check for hardcoreOverlap
                print("checking for hardcore Overlap")
                if s.quality(tolerance=0.3) > 0:
                    print("Found bad quality monomer insertion. Redoing last insertion...")
                    s.unwrap()
                    s.write_xyz('bad_insertion_' + str(insertion+2) + '.xyz')
                    s.wrap()
                    redoMonomerInsertion(s,n,insertion+2)
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
