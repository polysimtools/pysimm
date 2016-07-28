# ******************************************************************************
# pysimm.system module
# ******************************************************************************
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
import re
import sys
import json
from xml.etree import ElementTree as Et
from random import random
from StringIO import StringIO
from urllib2 import urlopen, HTTPError, URLError
from itertools import permutations
from math import sin, cos, sqrt, pi, acos, floor, ceil


try:
    from subprocess import call
except ImportError:
    call = None
try:
    import numpy as np
except ImportError:
    np = None

from pysimm import calc
from pysimm import error_print
from pysimm import warning_print
from pysimm import verbose_print
from pysimm import debug_print
from pysimm import PysimmError
from pysimm.calc import rotate_vector
from pysimm.utils import PysimmError, Item, ItemContainer


class Particle(Item):
    """pysimm.system.Particle

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        x: x coordinate
        y: y coordinate
        z: z coordinate
        charge: partial charge
        type: ParticleType object reference
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)

    def check(self, style='full'):
        if style == 'full':
            if self.x is None:
                error_print('particle %s has no x coordinate' % self.tag)
                return False
            if self.y is None:
                return False
            if self.z is None:
                return False
            if self.charge is None:
                return False
            if self.type is None or not self.type.check():
                return False
            return True
        else:
            error_print('style %s not supported yet' % style)
            return False

    def delete_bonding(self, s):
        """pysimm.system.Particle.delete_bonding

        Iterates through s.bonds, s.angles, s.dihedrals, and s.impropers and removes
        those which contain this Particle.

        Args:
            s: pysimm.system.System object from which bonding Item objects will be removed

        Returns:
            None
        """
        if self.bonds:
            for b in self.bonds:
                if b in s.bonds:
                    s.bonds.remove(b.tag)
        else:
            for b in s.bonds:
                if self is b.a or self is b.b:
                    s.bonds.remove(b.tag)

        if self.angles:
            for a in self.angles:
                if a in s.angles:
                    s.angles.remove(a.tag)
        else:
            for a in s.angles:
                if self is a.a or self is a.b or self is a.c:
                    s.angles.remove(a.tag)

        if self.dihedrals:
            for d in self.dihedrals:
                if d in s.dihedrals:
                    s.dihedrals.remove(d.tag)
        else:
            for d in s.dihedrals:
                if self is d.a or self is d.b or self is d.c or self is d.d:
                    s.dihedrals.remove(d.tag)

        if self.impropers:
            for i in self.impropers:
                if i in s.impropers:
                    s.impropers.remove(i.tag)
        else:
            for i in s.impropers:
                if self is i.a or self is i.b or self is i.c or self is i.d:
                    s.impropers.remove(i.tag)

    def __sub__(self, other):
        """pysimm.system.Particle.__sub__

        Implements subtraction between Particle objects to calculate distance.

        Args:
            other: pysimm.system.Particle object

        Returns:
            distance calculated by pysimm.calc.distance. This does not consider pbc
        """
        if isinstance(other, Particle):
            return calc.distance(self, other)
        else:
            return None

    def __rsub__(self, other):
        self.__sub__(other)


class ParticleType(Item):
    """pysimm.system.ParticleType

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        sigma: LJ sigma value (Angstrom)
        epsilon: LJ epsilon value (kcal/mol)
        elem: element abbreviation, i.e. 'H' for Hydrogen, 'Cl' for Chlorine
        name: force field particle type name
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)


class Bond(Item):
    """pysimm.system.Bond

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        a: pysimm.system.Particle object involved in bond
        b: pysimm.system.Particle object involved in bond
        type: BondType object reference
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)

    def distance(self):
        """pysimm.system.Bond.distance

        Calculates distance between Particle a and Particle b in this Bond object.
        Sets distance to dist attribute of Bond object. Does not consider pbc.

        Args:
            None

        Returns:
            Distance between Particle a and Particle b (not considering pbc)
        """
        if isinstance(self.a, Particle) and isinstance(self.b, Particle):
            self.dist = calc.distance(self.a, self.b)
            return self.dist
        else:
            return None


class BondType(Item):
    """pysimm.system.BondType

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        k: harmonic bond force constant (kcal/mol/A^2)
        r0: bond equilibrium distance (Angstrom)
        name: force field bond type name
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)
        if self.name:
            self.rname = ','.join(reversed(self.name.split(',')))


class Angle(Item):
    """pysimm.system.Angle

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        a: pysimm.system.Particle object involved in angle
        b: pysimm.system.Particle object involved in angle (middle particle)
        c: pysimm.system.Particle object involved in angle
        type: AngleType object reference
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)

    def angle(self, radians=False):
        """pysimm.system.Angle.angle

        Calculate angle.

        Args:
            radians: True to return value in radians (default: False)

        Returns:
            Angle between Particle a, b, and c
        """
        self.theta = calc.angle(self.a, self.b, self.c, radians)
        return self.theta


class AngleType(Item):
    """pysimm.system.AngleType

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        k: harmonic angle bend force constant (kcal/mol/radian^2)
        theta0: angle equilibrium value (degrees)
        name: force field angle type name
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)
        if self.name:
            self.rname = ','.join(reversed(self.name.split(',')))


class Dihedral(Item):
    """pysimm.system.Dihedral

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        a: pysimm.system.Particle object involved in dihedral
        b: pysimm.system.Particle object involved in dihedral (middle particle)
        c: pysimm.system.Particle object involved in dihedral (middle particle)
        d: pysimm.system.Particle object involved in dihedral
        type: DihedralType object reference
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)


class DihedralType(Item):
    """pysimm.system.DihedralType

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        k: dihedral energy barrier (kcal/mol)
        d: minimum (+1 or -1)
        n: multiplicity (integer >= 0)
        name: force field dihedral type name
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)
        if self.name:
            self.rname = ','.join(reversed(self.name.split(',')))


class Improper(Item):
    """pysimm.system.Improper

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        a: pysimm.system.Particle object involved in improper (middle particle)
        b: pysimm.system.Particle object involved in improper
        c: pysimm.system.Particle object involved in improper
        d: pysimm.system.Particle object involved in improper
        type: ImproperType object reference
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)


class ImproperType(Item):
    """pysimm.system.ImproperType

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        k: improper energy barrier (kcal/mol)
        x0: equilibrium value (degrees)
        name: force field improper type name
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)
        if self.name:
            self.rname = ','.join(reversed(self.name.split(',')))


class Dimension(Item):
    """pysimm.system.Dimension

    Objects inheriting from pysimm.utils.Item can contain arbitrary data.
    Keyword arguments are assigned as attributes.
    Attributes usually used are given below.

    Attributes:
        xlo: minimum value in x dimension
        xhi: maximum value in x dimension
        ylo: minimum value in y dimension
        yhi: maximum value in y dimension
        zlo: minimum value in z dimension
        zhi: maximum value in z dimension
        dx: distance in x dimension
        dy: distance in y dimension
        dz: distance in z dimension
    """
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)
        if (self.center and
                self.dx is not None and
                self.dy is not None and
                self.dz is not None):
            if self.center is True:
                self.center = [0., 0., 0.]
            self.xlo = -1*self.dx/2. + self.center[0]
            self.xhi = self.dx/2. + self.center[0]
            self.ylo = -1*self.dy/2. + self.center[1]
            self.yhi = self.dy/2. + self.center[1]
            self.zlo = -1*self.dz/2. + self.center[2]
            self.zhi = self.dz/2. + self.center[2]
        if self.xhi is not None and self.xlo is not None:
            self.dx = self.xhi - self.xlo
        if self.yhi is not None and self.ylo is not None:
            self.dy = self.yhi - self.ylo
        if self.zhi is not None and self.zlo is not None:
            self.dz = self.zhi - self.zlo

    def check(self):
        if ((self.xlo is not None and self.xhi is not None and
                self.ylo is not None and self.yhi is not None and
                self.zlo is not None and self.zhi is not None) and
                (self.dx is not None and self.dy is not None and
                 self.dz is not None)):
            return True
        elif self.center and self.dx and self.dy and self.dz:
            return True
        else:
            return False


class System(object):
    """pysimm.system.System

    Object representation of molecular system.
    Contains information required for molecular simulation.

    Attributes:
        dim: Dimension object reference
        particles: pysimm.utils.ItemContainer for Particle organization
        particle_types: pysimm.utils.ItemContainer for ParticleType organization
        bonds: pysimm.utils.ItemContainer for Bond organization
        bond_types: pysimm.utils.ItemContainer for BondType organization
        angles: pysimm.utils.ItemContainer for Angle organization
        angle_types: pysimm.utils.ItemContainer for AngleType organization
        dihedrals: pysimm.utils.ItemContainer for Dihedral organization
        dihedral_types: pysimm.utils.ItemContainer for DihedralType organization
        impropers: pysimm.utils.ItemContainer for Improper organization
        improper_types: pysimm.utils.ItemContainer for ImproperType organization
        molecules: pysimm.utils.ItemContainer for Molecule organization

    """
    def __init__(self, **kwargs):

        self.objectified = False

        self.name = kwargs.get('name') or 'pySIMM System Object'
        self.ff_class = kwargs.get('ff_class')
        self.ff_name = kwargs.get('ff_name')
        self.dim = Dimension(xlo=kwargs.get('xlo'), xhi=kwargs.get('xhi'),
                             ylo=kwargs.get('ylo'), yhi=kwargs.get('yhi'),
                             zlo=kwargs.get('zlo'), zhi=kwargs.get('zhi'),
                             dx=kwargs.get('dx'), dy=kwargs.get('dy'),
                             dz=kwargs.get('dz'), center=kwargs.get('center'))
        self.dim_check = self.dim.check()
        self.mass = kwargs.get('mass') or 0.0
        self.particle_types = kwargs.get('particle_types') or ItemContainer()
        self.bond_types = kwargs.get('bond_types') or ItemContainer()
        self.angle_types = kwargs.get('angle_types') or ItemContainer()
        self.dihedral_types = kwargs.get('dihedral_types') or ItemContainer()
        self.improper_types = kwargs.get('improper_types') or ItemContainer()
        self.molecule_types = kwargs.get('molecule_types') or ItemContainer()
        self.particles = kwargs.get('particles') or ItemContainer()
        self.bonds = kwargs.get('bonds') or ItemContainer()
        self.angles = kwargs.get('angles') or ItemContainer()
        self.dihedrals = kwargs.get('dihedrals') or ItemContainer()
        self.impropers = kwargs.get('impropers') or ItemContainer()
        self.molecules = kwargs.get('molecules') or ItemContainer()

        self.set_mass()
        self.set_volume()
        self.set_density()
        self.set_cog()

    def __getattr__(self, name):
        return None

    def copy(self, rotate_x=None, rotate_y=None, rotate_z=None,
             dx=0, dy=0, dz=0):
        """pysimm.system.System.copy

        Create duplicate System object. Default behavior does not modify particle positions.

        Args:
            rotate_x: rotate duplicate system around x axis by this value (radians)
            rotate_y: rotate duplicate system around y axis by this value (radians)
            rotate_z: rotate duplicate system around z axis by this value (radians)
            dx: translate duplicate system in x dimension by this value (Angstrom)
            dy: translate duplicate system in y dimension by this value (Angstrom)
            dz: translate duplicate system in z dimension by this value (Angstrom)
        """
        new = System()

        new.ff_class = self.ff_class
        new.ff_name = self.ff_name
        new.pair_style = self.pair_style

        new.dim = self.dim.copy()

        for _ in self.molecules:
            new.molecules.add(Molecule())

        for pt in self.particle_types:
            new.particle_types.add(pt.copy())

        for bt in self.bond_types:
            new.bond_types.add(bt.copy())

        for at in self.angle_types:
            new.angle_types.add(at.copy())

        for dt in self.dihedral_types:
            new.dihedral_types.add(dt.copy())

        for it in self.improper_types:
            new.improper_types.add(it.copy())

        for p in self.particles:
            new_p = p.copy()
            if p.type:
                new_p.type = new.particle_types[p.type.tag]
            new_p.molecule = new.molecules[p.molecule.tag]
            if rotate_x or rotate_y or rotate_z:
                new_p.x, new_p.y, new_p.z = rotate_vector(new_p.x, new_p.y, new_p.z,
                                                   rotate_x, rotate_y, rotate_z)
            new_p.x += dx
            new_p.y += dy
            new_p.z += dz

            new.particles.add(new_p)
            new_p.molecule.particles.add(new_p)
            new_p.bonds = ItemContainer()
            new_p.angles = ItemContainer()
            new_p.dihedrals = ItemContainer()
            new_p.impropers = ItemContainer()

        for b in self.bonds:
            new_b = b.copy()
            new_b.a=new.particles[b.a.tag]
            new_b.b=new.particles[b.b.tag]
            if b.type:
                new_b.type=new.bond_types[b.type.tag]
            new.bonds.add(new_b)
            new_b.a.molecule.bonds.add(new_b)
            new_b.a.bonds.add(new_b)
            new_b.b.bonds.add(new_b)

        for a in self.angles:
            new_a = Angle(a=new.particles[a.a.tag],
                          b=new.particles[a.b.tag],
                          c=new.particles[a.c.tag])
            if a.type:
                new_a.type=new.angle_types[a.type.tag]
            new.angles.add(new_a)
            new_a.a.molecule.angles.add(new_a)

        for d in self.dihedrals:
            new_d = Dihedral(a=new.particles[d.a.tag],
                             b=new.particles[d.b.tag],
                             c=new.particles[d.c.tag],
                             d=new.particles[d.d.tag])
            if d.type:
                new_d.type=new.dihedral_types[d.type.tag]
            new.dihedrals.add(new_d)
            new_d.a.molecule.dihedrals.add(new_d)

        for i in self.impropers:
            new_i = Improper(a=new.particles[i.a.tag],
                             b=new.particles[i.b.tag],
                             c=new.particles[i.c.tag],
                             d=new.particles[i.d.tag])
            if i.type:
                new_i.type=new.improper_types[i.type.tag]
            new.impropers.add(new_i)
            new_i.a.molecule.impropers.add(new_i)

        return new

    def add(self, other, **kwargs):
        """pysimm.system.System.add

        Add other System to this System. Optionally remove duplicate types (default behavior).

        Args:
            other: pysimm.system.System object to add to this System
            unique_types (optional): Remove duplicate types and reassign references to existing types (True)
            change_dim (optional): Update pysimm.system.Dimension object so that Particle objects do not exist
                                   outside of Dimension extremes (True)
            update_properties (optional): Update system-wide mass, volume, density, center of gravity, and velocity
                                          properties (True)
        """
        unique_types = kwargs.get('unique_types') if (kwargs.get('unique_types')
                                                      is not None) else True
        change_dim = kwargs.get('change_dim') if (kwargs.get('change_dim')
                                                  is not None) else True
        update_properties = kwargs.get('update_properties') if kwargs.get('update_properties') is not None else True

        if self.ff_class is not None and self.ff_class != other.ff_class:
            warning_print('warning: mixing forcefield classes is highly '
                          'unadvised. only continue if you know what you '
                          'are doing')

        for pt in other.particle_types:
            if unique_types:
                if pt.name not in [x.name for x in self.particle_types]:
                    del pt.tag
                    self.particle_types.add(pt)
            else:
                del pt.tag
                self.particle_types.add(pt)
        for bt in other.bond_types:
            if unique_types:
                if bt.name not in [x.name for x in self.bond_types]:
                    del bt.tag
                    self.bond_types.add(bt)
            else:
                del bt.tag
                self.bond_types.add(bt)
        for at in other.angle_types:
            if unique_types:
                if at.name not in [x.name for x in self.angle_types]:
                    del at.tag
                    self.angle_types.add(at)
            else:
                del at.tag
                self.angle_types.add(at)
        for dt in other.dihedral_types:
            if unique_types:
                if dt.name not in [x.name for x in self.dihedral_types]:
                    del dt.tag
                    self.dihedral_types.add(dt)
            else:
                del dt.tag
                self.dihedral_types.add(dt)
        for it in other.improper_types:
            if unique_types:
                if it.name not in [x.name for x in self.improper_types]:
                    del it.tag
                    self.improper_types.add(it)
            else:
                del it.tag
                self.improper_types.add(it)
        for p in other.particles:
            del p.tag
            if change_dim:
                self.dim.xhi = max(p.x, self.dim.xhi)
                self.dim.xlo = min(p.x, self.dim.xlo)
                self.dim.yhi = max(p.y, self.dim.yhi)
                self.dim.ylo = min(p.y, self.dim.ylo)
                self.dim.zhi = max(p.z, self.dim.zhi)
                self.dim.zlo = min(p.z, self.dim.zlo)
            if unique_types and p.type not in self.particle_types:
                pt = self.particle_types.get(p.type.name)
                if not pt or len(pt) > 1:
                    error_print('ParticleType error')
                else:
                    p.type = pt[0]
            self.particles.add(p)
        for b in other.bonds:
            del b.tag
            if unique_types and b.type not in self.bond_types:
                bt = self.bond_types.get(b.type.name)
                if not bt or len(bt) > 1:
                    error_print('BondType error')
                else:
                    b.type = bt[0]
            self.bonds.add(b)
        for a in other.angles:
            del a.tag
            if unique_types and a.type not in self.angle_types:
                at = self.angle_types.get(a.type.name)
                if not at or len(at) > 1:
                    error_print('AngleType error')
                else:
                    a.type = at[0]
            self.angles.add(a)
        for d in other.dihedrals:
            del d.tag
            if unique_types and d.type not in self.dihedral_types:
                dt = self.dihedral_types.get(d.type.name)
                if not dt:
                    error_print('DihedralType error')
                elif len(dt) > 1:
                    index = 0
                    x = 5
                    for i in range(len(dt)):
                        if dt[i].name.count('X') < x:
                            index = i
                            x = dt[i].name.count('X')
                    d.type = dt[index]
                else:
                    d.type = dt[0]
            self.dihedrals.add(d)
        for i in other.impropers:
            del i.tag
            if unique_types and i.type not in self.improper_types:
                it = self.improper_types.get(i.type.name)
                if not it or len(it) > 1:
                    error_print('ImproperType error')
                else:
                    i.type = it[0]
            self.impropers.add(i)

        for m in other.molecules:
            del m.tag
            self.molecules.add(m)

        if update_properties:
            self.set_mass()
            self.set_volume()
            self.set_density()
            self.set_cog()
            self.set_velocity()

    def distance(self, p1, p2):
        """pysimm.system.System.distance

        Calculate distance between two particles considering pbc.

        Args:
            p1: pysimm.system.Particle object
            p2: pysimm.system.Particle object

        Returns:
            distance between particles considering pbc
        """
        return calc.pbc_distance(self, p1, p2)

    def wrap(self):
        """pysimm.system.System.wrap

        Wrap Particle images into box defined by Dimension object.
        Ensure particles are contained within simulation box.

        Args:
            None

        Returns:
            None
        """
        self.dim.check()
        for p in self.particles:
            while p.x > self.dim.xhi:
                p.x -= self.dim.dx
            while p.x < self.dim.xlo:
                p.x += self.dim.dx
            while p.y > self.dim.yhi:
                p.y -= self.dim.dy
            while p.y < self.dim.ylo:
                p.y += self.dim.dy
            while p.z > self.dim.zhi:
                p.z -= self.dim.dz
            while p.z < self.dim.zlo:
                p.z += self.dim.dz

    def unwrap(self):
        """pysimm.system.System.unwrap()

        Unwraps Particle images such that no bonds cross box edges.

        Args:
            None

        Returns:
            None
        """
        self.dim.check()
        self.add_particle_bonding()
        self.add_particle_bonding()
        next_to_unwrap = []
        for p in self.particles:
            p.unwrapped = False
        for m in self.molecules:
            for p0 in m.particles:
                p0.unwrapped = True
                next_to_unwrap.append(p0)
                for p in next_to_unwrap:
                    for pb in p.bonded_to:
                        if pb.unwrapped:
                            continue
                        next_to_unwrap.append(pb)
                        pb.unwrapped = True
                        dx = p.x - pb.x
                        while abs(dx) > self.dim.dx / 2:
                            if dx > 0:
                                pb.x += self.dim.dx
                            else:
                                pb.x -= self.dim.dx
                            dx = p.x - pb.x
                        dy = p.y - pb.y
                        while abs(dy) > self.dim.dy / 2:
                            if dy > 0:
                                pb.y += self.dim.dy
                            else:
                                pb.y -= self.dim.dy
                            dy = p.y - pb.y
                        dz = p.z - pb.z
                        while abs(dz) > self.dim.dz / 2:
                            if dz > 0:
                                pb.z += self.dim.dz
                            else:
                                pb.z -= self.dim.dz
                            dz = p.z - pb.z

        for b in self.bonds:
            if b.distance() > 5:
                print('unwrap probably failed')
                return False

        return True

    def quality(self, tolerance=0.1):
        """pysimm.system.System.quality

        Attemps to assess quality of System based on bond lengths in unwrapped system.

        Args:
            tolerance: fractional value of equilibrium bond length that is acceptable

        Returns:
            number of bonds in system outside tolerance
        """
        self.unwrap()
        bad_bonds = 0
        for b in self.bonds:
            if b.distance() > b.type.r0*(1+tolerance) or b.distance() < b.type.r0*(1-tolerance):
                bad_bonds += 1
        verbose_print('%s of %s bonds found to be outside of tolerance' % (bad_bonds, self.bonds.count))
        self.wrap()

    def shift_to_origin(self):
        """pysimm.system.System.shift_to_origin

        Shifts simulation box to begin at origin. i.e. xlo=ylo=zlo=0

        Args:
            None

        Returns:
            None
        """
        for p in self.particles:
            p.x -= self.dim.xlo
            p.y -= self.dim.ylo
            p.z -= self.dim.zlo
        self.dim.xhi -= self.dim.xlo
        self.dim.yhi -= self.dim.ylo
        self.dim.zhi -= self.dim.zlo
        self.dim.xlo -= self.dim.xlo
        self.dim.ylo -= self.dim.ylo
        self.dim.zlo -= self.dim.zlo

    def set_charge(self):
        """pysimm.system.System.set_charge

        Sets total charge of all Particle objects in System.particles

        Args:
            None

        Returns:
            None
        """
        self.charge = 0
        for p in self.particles:
            self.charge += p.charge

    def zero_charge(self):
        """pysimm.system.System.zero_charge

        Enforces total System charge to be 0.0 by subtracting excess charge from last particle

        Args:
            None:

        Returns:
            None
        """
        charge = 0.
        for p in self.particles:
            charge += p.charge
        if charge != 0:
            p.charge -= charge
        self.set_charge()

    def check_items(self):
        if len(self.particles) != self.particles.count:
            raise PysimmError('particles missing')
        if len(self.bonds) != self.bonds.count:
            raise PysimmError('bonds missing')
        if len(self.angles) != self.angles.count:
            raise PysimmError('angles missing')
        if len(self.dihedrals) != self.dihedrals.count:
            raise PysimmError('dihedrals missing')
        if len(self.impropers) != self.impropers.count:
            raise PysimmError('impropers missing')
        if len(self.molecules) != self.molecules.count:
            raise PysimmError('molecules missing')

    def update_particle_types_from_forcefield(self, f):
        """pysimm.system.System.update_types_from_forcefield

        Updates ParticleType data from Forcefield object f based on ParticleType.name

        Args:
            f: pysimm.forcefield.Forcefield object reference

        Returns:
            None
        """
        for pt in self.particle_types:
            name_ = pt.name.split('@')[-1]
            linker = False
            if pt.name.find('@') >= 0:
                linker = pt.name.split('@')[0]
            pt_ = f.particle_types.get(name_)
            if pt_:
                new = pt_[0].copy()
                new.tag = pt.tag
                if linker:
                    new.name = '%s@%s' % (linker, new.name)
                self.particle_types.remove(pt.tag)
                self.particle_types.add(new)

    def make_linker_types(self):
        """pysimm.system.System.make_linker_types

        Identifies linker particles and creates duplicate ParticleType objects with new names.
        Identification is performed by Particle.linker attribute.
        New ParticleType name is prepended with [H or T]L@ to designate head or tail linker

        Args:
            None

        Returns:
            None
        """
        for p in self.particles:
            if p.linker == 'head':
                head_linker = self.particle_types.get('HL@%s' % p.type.name)
                if head_linker:
                    p.type = head_linker[0]
                else:
                    p.type = p.type.copy()
                    p.type.name = 'HL@%s' % p.type.name
                    self.particle_types.add(p.type)
            elif p.linker == 'tail':
                tail_linker = self.particle_types.get('TL@%s' % p.type.name)
                if tail_linker:
                    p.type = tail_linker[0]
                else:
                    p.type = p.type.copy()
                    p.type.name = 'TL@%s' % p.type.name
                    self.particle_types.add(p.type)
            elif p.linker:
                linker = self.particle_types.get('L@%s' % p.type.name)
                if linker:
                    p.type = linker[0]
                else:
                    p.type = p.type.copy()
                    p.type.name = 'L@%s' % p.type.name
                    self.particle_types.add(p.type)

    def remove_linker_types(self):
        """pysimm.system.System.remove_linker_types

        Reassigns Particle.type references to original ParticleType objects without linker prepend

        Args:
            None

        Returns:
            None
        """
        for p in self.particles:
            if p.type.name.find('@') >= 0:
                pt = self.particle_types.get(p.type.name.split('@')[-1])
                if pt:
                    p.type = pt[0]
                else:
                    print('cannot find regular type for linker %s'
                          % p.type.name)

    def read_lammpstrj(self, trj, frame=1):
        """pysimm.system.System.read_lammpstrj

        Updates particle positions and box size from LAMMPS trajectory file at given frame

        Args:
            trj: LAMMPS trajectory file
            frame: sequential frame number (not LAMMPS timestep) default=1

        Returns:
            None
        """
        t_frame = 0
        nparticles = 0
        updated = 0
        with open(trj) as f:
            line = f.readline()
            while line:
                if len(line.split()) > 1 and line.split()[1] == 'TIMESTEP':
                    t_frame += 1
                elif len(line.split()) > 1 and line.split()[1] == 'NUMBER':
                    nparticles = int(f.readline())
                elif len(line.split()) > 1 and line.split()[1] == 'BOX':
                    self.dim.xlo, self.dim.xhi = map(float,
                                                     f.readline().split())
                    self.dim.ylo, self.dim.yhi = map(float,
                                                     f.readline().split())
                    self.dim.zlo, self.dim.zhi = map(float,
                                                     f.readline().split())
                elif (len(line.split()) > 1 and line.split()[1] == 'ATOMS' and
                      t_frame == frame):
                    for i in range(nparticles):
                        line = f.readline().split()
                        if len(line) == 4:
                            id_, x, y, z = map(float, line)
                        elif len(line) == 5:
                            id_, type_, x, y, z = map(float, line)
                        elif len(line) == 7:
                            id_, type_, x, y, z, ix, iy, iz = map(float, line)
                        else:
                            error_print('cannot understand lammpstrj formatting; exiting')
                            return
                        id_ = int(id_)
                        if self.particles[id_]:
                            updated += 1
                            self.particles[id_].x = x
                            self.particles[id_].y = y
                            self.particles[id_].z = z
                line = f.readline()

        verbose_print('updated particle positions for %s of %s particles from trajectory' % (updated, nparticles))

    def read_xyz(self, xyz, frame=1, quiet=False):
        """pysimm.system.System.read_xyz

        Updates particle positions and box size from xyz file at given frame

        Args:
            xyz: xyz trajectory file
            frame: sequential frame number default=1
            quiet: True to print status default=False

        Returns:
            None
        """
        if not quiet:
            verbose_print('reading particle positions from %s' % xyz)
            warning_print('particles are assumed to be in order in xyz file')
        t_frame = 0
        with open(xyz) as f:
            line = f.readline()
            while line:
                t_frame += 1
                assert int(line.split()[0]) == self.particles.count
                line = f.readline()
                for n in range(1, self.particles.count + 1):
                    p = self.particles[n]
                    if t_frame == 1:
                        if not p.type.elem and p.type.name:
                            if p.type.name[0].lower() != 'l':
                                p.type.elem = p.type.name[0].upper()
                            else:
                                p.type.elem = p.type.name[1].upper()
                    line = f.readline()
                    if t_frame == frame:
                        x, y, z = map(float, line.split()[1:])
                        p.x = x
                        p.y = y
                        p.z = z
                if t_frame == frame:
                    print('read %s particle positions from %s'
                          % (self.particles.count, xyz))
                line = f.readline()

    def update_types(self, ptypes, btypes, atypes, dtypes, itypes):
        """pysimm.system.System.update_types

        Updates type objects from a given list of types.

        Args:
            ptypes: list of ParticleType objects from which to update
            btypes: list of BondType objects from which to update
            atypes: list of AngleType objects from which to update
            dtypes: list of DihedralType objects from which to update
            itypes: list of ImproperType objects from which to update
        """
        if ptypes is not None:
            for p in self.particles:
                pt = self.particle_types.get(p.type.name, first=True)
                if pt:
                    p.type = pt[0]
            self.particle_types.remove('all')
            for pt in ptypes:
                self.particle_types.add(pt)

        if btypes is not None:
            for b in self.bonds:
                bt = self.bond_types.get(b.type.name, first=True)
                if bt:
                    b.type = bt[0]
            self.bond_types.remove('all')
            for bt in btypes:
                self.bond_types.add(bt)

        if atypes is not None:
            for a in self.angles:
                at = self.angle_types.get(a.type.name, first=True)
                if at:
                    a.type = at[0]
            self.angle_types.remove('all')
            for at in atypes:
                self.angle_types.add(at)

        if dtypes is not None:
            for d in self.dihedrals:
                dt = self.dihedral_types.get(d.type.name, first=True)
                if dt:
                    d.type = dt[0]
            self.dihedral_types.remove('all')
            for dt in dtypes:
                self.dihedral_types.add(dt)

        if itypes is not None:
            for i in self.impropers:
                it = self.improper_types.get(i.type.name, first=True)
                if it:
                    i.type = it[0]
            self.improper_types.remove('all')
            for it in itypes:
                self.improper_types.add(it)

    def read_type_names(self, types_file):
        """pysimm.system.System.read_type_names

        Update ParticleType names from file.

        Args:
            types_file: type dictionary file name

        Returns:
            None
        """
        ptypes = dict()
        btypes = dict()
        atypes = dict()
        dtypes = dict()
        itypes = dict()

        if os.path.isfile(types_file):
            f = file(types_file)
        elif isinstance(types_file, basestring):
            f = StringIO(types_file)
        for line in f:
            line = line.split()
            if line and line[0].lower() == 'atom':
                for i in range(self.particle_types.count):
                    line = f.next().split()
                    ptypes[int(line[0])] = line[1]
            elif line and line[0].lower() == 'bond':
                for i in range(self.bond_types.count):
                    line = f.next().split()
                    btypes[int(line[0])] = line[1]
            elif line and line[0].lower() == 'angle':
                for i in range(self.angle_types.count):
                    line = f.next().split()
                    atypes[int(line[0])] = line[1]
            elif line and line[0].lower() == 'dihedral':
                for i in range(self.dihedral_types.count):
                    line = f.next().split()
                    dtypes[int(line[0])] = line[1]
            elif line and line[0].lower() == 'improper':
                for i in range(self.improper_types.count):
                    line = f.next().split()
                    itypes[int(line[0])] = line[1]

        for t in self.particle_types:
            t.name = ptypes[t.tag]
            if t.name[0] == 'L':
                if t.name[1].upper() in ['H', 'C', 'N', 'O']:
                    t.elem = t.name[1].upper()
            else:
                if t.name[0].upper() in ['H', 'C', 'N', 'O']:
                    t.elem = t.name[0].upper()
        for t in self.bond_types:
            t.name = btypes[t.tag]
            t.rname = ','.join(reversed(t.name.split(',')))
        for t in self.angle_types:
            t.name = atypes[t.tag]
            t.rname = ','.join(reversed(t.name.split(',')))
        for t in self.dihedral_types:
            t.name = dtypes[t.tag]
            t.rname = ','.join(reversed(t.name.split(',')))
        for t in self.improper_types:
            t.name = itypes[t.tag]
            t.rname = ','.join(reversed(t.name.split(',')))

    def remove_spare_bonding(self, update_tags=True):
        """pysimm.system.System.remove_spare_bonding

        Removes bonds, angles, dihedrals and impropers that reference particles not in System.particles

        Args:
            update_tags: True to update all tags after removal of bonding items default=True
        """
        for b in self.bonds:
            if b.a not in self.particles or b.b not in self.particles:
                self.bonds.remove(b.tag, update=False)

        for a in self.angles:
            if (a.a not in self.particles or a.b not in self.particles or
                    a.c not in self.particles):
                self.angles.remove(a.tag, update=False)

        for d in self.dihedrals:
            if (d.a not in self.particles or d.b not in self.particles or
                    d.c not in self.particles or d.d not in self.particles):
                self.dihedrals.remove(d.tag, update=False)

        for i in self.impropers:
            if (i.a not in self.particles or i.b not in self.particles or
                    i.c not in self.particles or i.d not in self.particles):
                self.impropers.remove(i.tag, update=False)

        if update_tags:
            self.update_tags()

    def update_tags(self):
        """pysimm.system.System.update_tags

        Update Item tags in ItemContainer objects to preserve continuous tags

         Args:
             None

         Returns:
             None
        """
        particles = self.particles.get('all')
        self.particles.remove('all')
        for p in particles:
            del p.tag
            self.particles.add(p)

        ptypes = self.particle_types.get('all')
        self.particle_types.remove('all')
        for pt in ptypes:
            del pt.tag
            self.particle_types.add(pt)

        bonds = self.bonds.get('all')
        self.bonds.remove('all')
        for b in bonds:
            del b.tag
            self.bonds.add(b)

        btypes = self.bond_types.get('all')
        self.bond_types.remove('all')
        for bt in btypes:
            del bt.tag
            self.bond_types.add(bt)

        angles = self.angles.get('all')
        self.angles.remove('all')
        for a in angles:
            del a.tag
            self.angles.add(a)

        atypes = self.angle_types.get('all')
        self.angle_types.remove('all')
        for at in atypes:
            del at.tag
            self.angle_types.add(at)

        dihedrals = self.dihedrals.get('all')
        self.dihedrals.remove('all')
        for d in dihedrals:
            del d.tag
            self.dihedrals.add(d)

        dtypes = self.dihedral_types.get('all')
        self.dihedral_types.remove('all')
        for dt in dtypes:
            del dt.tag
            self.dihedral_types.add(dt)

        impropers = self.impropers.get('all')
        self.impropers.remove('all')
        for i in impropers:
            del i.tag
            self.impropers.add(i)

        itypes = self.improper_types.get('all')
        self.improper_types.remove('all')
        for it in itypes:
            del it.tag
            self.improper_types.add(it)

    def set_references(self):
        """pysimm.system.System.set_references

        Set object references when System information read from text file.
        For example, if bond type value 2 is read from file, set Bond.type to bond_types[2]

        Args:
            None

        Returns:
            None
        """
        for p in self.particles:
            if isinstance(p.type, int) and self.particle_types[p.type]:
                p.type = self.particle_types[p.type]
            elif isinstance(p.type, int) and not self.particle_types[p.type]:
                error_print('error: Cannot find type with tag %s in system '
                            'particles types' % p.type)

        for b in self.bonds:
            if isinstance(b.type, int) and self.bond_types[b.type]:
                b.type = self.bond_types[b.type]
            elif isinstance(b.type, int) and not self.bond_types[b.type]:
                error_print('error: Cannot find type with tag %s in system '
                            'bond types' % b.type)

        for a in self.angles:
            if isinstance(a.type, int) and self.angle_types[a.type]:
                a.type = self.angle_types[a.type]
            elif isinstance(b.type, int) and not self.angle_types[a.type]:
                error_print('error: Cannot find type with tag %s in system '
                            'angle types' % a.type)

        for d in self.dihedrals:
            if isinstance(d.type, int) and self.dihedral_types[d.type]:
                d.type = self.dihedral_types[d.type]
            elif isinstance(d.type, int) and not self.dihedral_types[d.type]:
                error_print('error: Cannot find type with tag %s in system '
                            'angle types' % d.type)

        for i in self.impropers:
            if isinstance(i.type, int) and self.improper_types[i.type]:
                i.type = self.improper_types[i.type]
            elif isinstance(i.type, int) and not self.improper_types[i.type]:
                error_print('error: Cannot find type with tag %s in system '
                            'angle types' % i.type)

    def objectify(self):
        """pysimm.system.System.objectify

        Set references for Bond, Angle, Dihedral, Improper objects.
        For example, if read from file that bond #1 is between particle 1 and 2 set Bond.a to particles[1], etc.

        Args:
            None

        Returns:
            None
        """
        if self.objectified:
            return 'already objectified'
        self.set_references()
        for p in self.particles:
            if not isinstance(p.molecule, Molecule):
                if not self.molecules[p.molecule]:
                    m = Molecule()
                    m.tag = p.molecule
                    self.molecules.add(m)
                p.molecule = self.molecules[p.molecule]
                self.molecules[p.molecule.tag].particles.add(p)
            p.bonds = ItemContainer()
            p.angles = ItemContainer()
            p.dihedrals = ItemContainer()
            p.impropers = ItemContainer()
        for b in self.bonds:
            if type(b.a) == int:
                b.a = self.particles[b.a]
                b.b = self.particles[b.b]
                b.a.bonds.add(b)
                b.b.bonds.add(b)
                if b.a.molecule:
                    b.a.molecule.bonds.add(b)
        for a in self.angles:
            if type(a.a) == int:
                a.a = self.particles[a.a]
                a.b = self.particles[a.b]
                a.c = self.particles[a.c]
                if a.a.molecule:
                    a.a.molecule.angles.add(a)
        for d in self.dihedrals:
            if type(d.a) == int:
                d.a = self.particles[d.a]
                d.b = self.particles[d.b]
                d.c = self.particles[d.c]
                d.d = self.particles[d.d]
                if d.a.molecule:
                    d.a.molecule.dihedrals.add(d)
        for i in self.impropers:
            if type(i.a) == int:
                i.a = self.particles[i.a]
                i.b = self.particles[i.b]
                i.c = self.particles[i.c]
                i.d = self.particles[i.d]
                if i.a.molecule:
                    i.a.molecule.impropers.add(i)
        self.objectified = True

    def add_particle_bonding(self):
        """pysimm.system.System.add_particle_bonding

        Update Particle objects such that Particle.bonded_to contains other Particle objects invloved in bonding

        Args:
            None

        Returns:
            None
        """
        for p in self.particles:
            p.bonded_to = ItemContainer()
            p.bonds = ItemContainer()
        for b in self.bonds:
            b.a.bonded_to.add(b.b)
            b.a.bonds.add(b)
            b.b.bonded_to.add(b.a)
            b.b.bonds.add(b)

    def set_excluded_particles(self, bonds=True, angles=True, dihedrals=True):
        """pysimm.system.System.set_excluded_particles

        Updates Particle object such that Particle.excluded_particles contains other Particle objects involved in
        1-2, 1-3, and/or 1-4 interactions

        Args:
            bonds: exclude particles involved in 1-2 interactions
            angles: exclude particles involved in 1-3 interactions
            dihedrals: exclude particles involved in 1-4 interactions

        """
        for p in self.particles:
            p.excluded_particles = ItemContainer()

        if bonds:
            for b in self.bonds:
                if b.a.tag < b.b.tag:
                    b.a.excluded_particles.add(b.b)
                else:
                    b.b.excluded_particles.add(b.a)

        if angles:
            for a in self.angles:
                if a.a.tag < a.b.tag:
                    a.a.excluded_particles.add(a.b)
                if a.a.tag < a.c.tag:
                    a.a.excluded_particles.add(a.c)
                if a.b.tag < a.a.tag:
                    a.b.excluded_particles.add(a.a)
                if a.b.tag < a.c.tag:
                    a.b.excluded_particles.add(a.c)
                if a.c.tag < a.a.tag:
                    a.c.excluded_particles.add(a.a)
                if a.c.tag < a.b.tag:
                    a.c.excluded_particles.add(a.b)

        if dihedrals:
            for d in self.dihedrals:
                if d.a.tag < d.b.tag:
                    d.a.excluded_particles.add(d.b)
                if d.a.tag < d.c.tag:
                    d.a.excluded_particles.add(d.c)
                if d.a.tag < d.d.tag:
                    d.a.excluded_particles.add(d.d)
                if d.b.tag < d.a.tag:
                    d.b.excluded_particles.add(d.a)
                if d.b.tag < d.c.tag:
                    d.b.excluded_particles.add(d.c)
                if d.b.tag < d.d.tag:
                    d.b.excluded_particles.add(d.d)
                if d.c.tag < d.a.tag:
                    d.c.excluded_particles.add(d.a)
                if d.c.tag < d.b.tag:
                    d.c.excluded_particles.add(d.b)
                if d.c.tag < d.d.tag:
                    d.c.excluded_particles.add(d.d)
                if d.d.tag < d.a.tag:
                    d.d.excluded_particles.add(d.a)
                if d.d.tag < d.b.tag:
                    d.d.excluded_particles.add(d.b)
                if d.d.tag < d.c.tag:
                    d.d.excluded_particles.add(d.c)

    def set_atomic_numbers(self):
        """pysimm.system.System.set_atomic_numbers

        Updates ParticleType objects with atomic number based on ParticleType.elem

        Args:
            None

        Returns:
            None
        """
        for pt in self.particle_types:
            if pt.elem == 'H':
                pt.atomic_number = 1
            elif pt.elem == 'He':
                pt.atomic_number = 2
            elif pt.elem == 'Li':
                pt.atomic_number = 3
            elif pt.elem == 'Be':
                pt.atomic_number = 4
            elif pt.elem == 'B':
                pt.atomic_number = 5
            elif pt.elem == 'C':
                pt.atomic_number = 6
            elif pt.elem == 'N':
                pt.atomic_number = 7
            elif pt.elem == 'O':
                pt.atomic_number = 8
            elif pt.elem == 'F':
                pt.atomic_number = 9
            elif pt.elem == 'Ne':
                pt.atomic_number = 10
            elif pt.elem == 'Na':
                pt.atomic_number = 11
            elif pt.elem == 'Mg':
                pt.atomic_number = 12
            elif pt.elem == 'Al':
                pt.atomic_number = 13
            elif pt.elem == 'Si':
                pt.atomic_number = 14
            elif pt.elem == 'P':
                pt.atomic_number = 15
            elif pt.elem == 'S':
                pt.atomic_number = 16
            elif pt.elem == 'Cl':
                pt.atomic_number = 17
            elif pt.elem == 'Ar':
                pt.atomic_number = 18
            elif pt.elem == 'K':
                pt.atomic_number = 19
            elif pt.elem == 'Ca':
                pt.atomic_number = 20
            elif pt.elem == 'Br':
                pt.atomic_number = 35
            elif pt.elem == 'I':
                pt.atomic_number = 53

    def add_particle_bonded_to(self, p, p0, f=None, sep=1.5):
        """pysimm.system.System.add_particle_bonded_to

        Add new Particle to System bonded to p0 and automatically update new forcefield types

        Args:
            p: new Particle object to be added
            p0: original Particle object in System to which p will be bonded
            f: pysimm.forcefield.Forcefield object from which new force field types will be retrieved

        Returns:
            new Particle being added to system for convenient reference
        """
        if p.x is None or p.y is None or p.z is None:
            phi = random() * 2 * pi
            theta = acos(random() * 2 - 1)
            p.x = p0.x + sep * cos(theta) * sin(phi)
            p.y = p0.y + sep * sin(theta) * sin(phi)
            p.z = p0.z + sep * cos(phi)
        if p.charge is None:
            p.charge = 0
        if p.molecule is None:
            p.molecule = p0.molecule
        self.add_particle(p)
        self.add_bond(p0, p, f)
        if not p0.bonded_to:
            self.add_particle_bonding()
        for p_ in p0.bonded_to:
            if p_ is not p:
                self.add_angle(p_, p0, p, f)
                for p_b in p_.bonded_to:
                    if p_b is not p0:
                        self.add_dihedral(p_b, p_, p0, p, f)
        return p

    def add_particle(self, p):
        self.particles.add(p)

    def rotate(self, around=None, theta_x=0, theta_y=0, theta_z=0, rot_matrix=None):
        """pysimm.system.System.rotate

        *** REQUIRES NUMPY ***

        Rotates entire system around given Particle by user defined angles

        Args:
            around: Particle around which System will be rotated default=None
            theta_x: angle around which system will be rotated on x axis
            theta_y: angle around which system will be rotated on y axis
            theta_z: angle around which system will be rotated on z axis
            rot_matrix: rotation matrix to use for rotation

        Returns:
            None
        """
        if not np:
            raise PysimmError('pysimm.system.System.rotate function requires numpy')
        theta_x = random() * 2 * pi if theta_x is 'random' else theta_x
        theta_y = random() * 2 * pi if theta_y is 'random' else theta_y
        theta_z = random() * 2 * pi if theta_z is 'random' else theta_z
        if around is None:
            around = []
            self.set_cog()
            around.append(self.cog[0])
            around.append(self.cog[1])
            around.append(self.cog[2])
        elif isinstance(around, Particle):
            around = [around.x, around.y, around.z]
        if (isinstance(around, list) and len(around) == 3 and
                len(set([isinstance(x, float) for x in around])) == 1 and isinstance(around[0], float)):
            for p in self.particles:
                p.x -= around[0]
                p.y -= around[1]
                p.z -= around[2]
                if rot_matrix is not None:
                    p.x, p.y, p.z = [x[0] for x in (rot_matrix*np.matrix([[p.x], [p.y], [p.z]])).tolist()]
                else:
                    p.x, p.y, p.z = rotate_vector(p.x, p.y, p.z, theta_x, theta_y, theta_z)
                p.x += around[0]
                p.y += around[1]
                p.z += around[2]

    def make_new_bonds(self, p1=None, p2=None, f=None, angles=True, dihedrals=True, impropers=True):
        """pysimm.system.System.make_new_bonds

        Makes new bond between two particles and updates new force field types

        Args:
            p1: pysimm.system.Particle object involved in new bond
            p2: pysimm.system.Particle object involved in new bond
            f: pysimm.forcefield.Forcefield object from which new force field types will be retrieved
            angles: True to update new angles default=True
            dihedrals: True to update new dihedrals default=True
            impropers: True to update new impropers default=True

        Returns:
            None
        """
        self.add_particle_bonding()
        self.add_bond(p1, p2, f)
        if angles or dihedrals or impropers:
            for p in p1.bonded_to:
                if angles:
                    self.add_angle(p, p1, p2, f)
                if dihedrals:
                    for pb in p.bonded_to:
                        if pb is not p1:
                            self.add_dihedral(pb, p, p1, p2, f)
            for p in p2.bonded_to:
                if angles:
                    self.add_angle(p1, p2, p, f)
                if dihedrals:
                    for pb in p.bonded_to:
                        if pb is not p2:
                            self.add_dihedral(p1, p2, p, pb, f)
            if dihedrals:
                for pb1 in p1.bonded_to:
                    for pb2 in p2.bonded_to:
                        self.add_dihedral(pb1, p1, p2, pb2, f)

        p1.bonded_to.add(p2)
        p2.bonded_to.add(p1)

        if impropers:
            if self.ff_class == '2':
                for perm in permutations(p1.bonded_to, 3):
                    unique = True
                    for i in self.impropers:
                        if i.a is not p1:
                            continue
                        if set([i.b, i.c, i.d]) == set([perm[0], perm[1],
                                                        perm[2]]):
                            unique = False
                            break
                    if unique:
                        self.add_improper(p1, perm[0], perm[1], perm[2], f)
                for perm in permutations(p2.bonded_to, 3):
                    unique = True
                    for i in self.impropers:
                        if i.a is not p2:
                            continue
                        if set([i.b, i.c, i.d]) == set([perm[0], perm[1],
                                                        perm[2]]):
                            unique = False
                            break
                    if unique:
                        self.add_improper(p2, perm[0], perm[1], perm[2], f)

    def add_bond(self, a=None, b=None, f=None):
        """pysimm.system.System.add_bond

        Add Bond to system between two particles

        Args:
            a: pysimm.system.Particle involved in new Bond
            b: pysimm.system.Particle involved in new Bond
            f: pysimm.forcefield.Forcefield object from which new force field type will be retrieved

        Returns:
            None
        """
        a_name = a.type.eq_bond or a.type.name
        b_name = b.type.eq_bond or b.type.name
        btype = self.bond_types.get('%s,%s' % (a_name, b_name))
        if not btype and f:
            btype = f.bond_types.get('%s,%s' % (a_name, b_name))
            if btype:
                bt = btype[0].copy()
                self.bond_types.add(bt)
            btype = self.bond_types.get('%s,%s' % (a_name, b_name))
        if btype:
            new_b = Bond(type=btype[0], a=a, b=b)
            self.bonds.add(new_b)
            if a.bonded_to is None or b.bonded_to is None:
                self.add_particle_bonding()
            if a.bonded_to and b not in a.bonded_to:
                a.bonded_to.add(b)
            if b.bonded_to and a not in b.bonded_to:
                b.bonded_to.add(a)
        else:
            error_print('error: system does not contain bond type named %s,%s '
                        'or could not find type in forcefield supplied'
                        % (a_name, b_name))
            return

    def add_angle(self, a=None, b=None, c=None, f=None):
        """pysimm.system.System.add_angle

        Add Angle to system between three particles

        Args:
            a: pysimm.system.Particle involved in new Angle
            b: pysimm.system.Particle involved in new Angle (middle particle)
            c: pysimm.system.Particle involved in new Angle
            f: pysimm.forcefield.Forcefield object from which new force field type will be retrieved

        Returns:
            None
        """
        a_name = a.type.eq_angle or a.type.name
        b_name = b.type.eq_angle or b.type.name
        c_name = c.type.eq_angle or c.type.name
        atype = self.angle_types.get('%s,%s,%s'
                                     % (a_name, b_name, c_name))
        if not atype and f:
            atype = f.angle_types.get('%s,%s,%s'
                                      % (a_name, b_name, c_name))
            if atype:
                at = atype[0].copy()
                self.angle_types.add(at)
            atype = self.angle_types.get('%s,%s,%s'
                                         % (a_name, b_name,
                                            c_name))
        if atype:
            self.angles.add(Angle(type=atype[0], a=a, b=b, c=c))
        else:
            error_print('error: system does not contain angle type named '
                        '%s,%s,%s or could not find type in forcefield supplied'
                        % (a_name, b_name, c_name))
            return

    def add_dihedral(self, a=None, b=None, c=None, d=None, f=None):
        """pysimm.system.System.add_dihedral

        Add Dihedral to system between four particles

        Args:
            a: pysimm.system.Particle involved in new Dihedral
            b: pysimm.system.Particle involved in new Dihedral (middle particle)
            c: pysimm.system.Particle involved in new Dihedral (middle particle)
            d: pysimm.system.Particle involved in new Dihedral
            f: pysimm.forcefield.Forcefield object from which new force field type will be retrieved

        Returns:
            None
        """
        a_name = a.type.eq_dihedral or a.type.name
        b_name = b.type.eq_dihedral or b.type.name
        c_name = c.type.eq_dihedral or c.type.name
        d_name = d.type.eq_dihedral or d.type.name
        dtype = self.dihedral_types.get('%s,%s,%s,%s'
                                        % (a_name, b_name,
                                           c_name, d_name))
        if not dtype and f:
            dtype = f.dihedral_types.get('%s,%s,%s,%s'
                                         % (a_name, b_name,
                                            c_name, d_name))
            if dtype:
                dt = dtype[0].copy()
                self.dihedral_types.add(dt)
            dtype = self.dihedral_types.get('%s,%s,%s,%s'
                                            % (a_name, b_name,
                                               c_name, d_name))
        if dtype:
            self.dihedrals.add(Dihedral(type=dtype[0], a=a, b=b, c=c, d=d))
        else:
            error_print('error: system does not contain dihedral type named '
                        '%s,%s,%s,%s or could not find type in forcefield '
                        'supplied' % (a_name, b_name,
                                      c_name, d_name))
            error_print('tags: %s %s %s %s' % (a.tag, b.tag, c.tag, d.tag))
            return

    def add_improper(self, a=None, b=None, c=None, d=None, f=None):
        """pysimm.system.System.add_improper

        Add Improper to system between four particles

        Args:
            a: pysimm.system.Particle involved in new Improper (middle particle)
            b: pysimm.system.Particle involved in new Improper
            c: pysimm.system.Particle involved in new Improper
            d: pysimm.system.Particle involved in new Improper
            f: pysimm.forcefield.Forcefield object from which new force field type will be retrieved

        Returns:
            None
        """
        a_name = a.type.eq_improper or a.type.name
        b_name = b.type.eq_improper or b.type.name
        c_name = c.type.eq_improper or c.type.name
        d_name = d.type.eq_improper or d.type.name
        if self.ff_class == '2' or self.improper_style == 'class2':
            itype = self.improper_types.get('%s,%s,%s,%s'
                                            % (b_name, a_name,
                                               c_name, d_name),
                                            improper_type=True)
        else:
            itype = self.improper_types.get('%s,%s,%s,%s'
                                            % (a_name, b_name,
                                               c_name, d_name),
                                            improper_type=True)
        if not itype and f:
            if f.ff_class == '2':
                itype = f.improper_types.get('%s,%s,%s,%s'
                                             % (b_name, a_name,
                                                c_name, d_name),
                                             improper_type=True)
            else:
                itype = f.improper_types.get('%s,%s,%s,%s'
                                             % (a_name, b_name,
                                                c_name, d_name),
                                             improper_type=True)
            if itype:
                it = itype[0].copy()
                self.improper_types.add(it)
            if f.ff_class == '2':
                itype = self.improper_types.get('%s,%s,%s,%s'
                                                % (b_name, a_name,
                                                   c_name, d_name),
                                                improper_type=True)
            else:
                itype = self.improper_types.get('%s,%s,%s,%s'
                                                % (a_name, b_name,
                                                   c_name, d_name),
                                                improper_type=True)
        if itype:
            self.impropers.add(Improper(type=itype[0], a=a, b=b, c=c, d=d))
        else:
            return

    def check_forcefield(self):
        if not self.objectified:
            self.objectify()
        for p in self.particles:
            if not p.bond_elements:
                p.bond_elements = [x.a.type.elem if p is x.b else
                                   x.b.type.elem for x in p.bonds]
                p.nbonds = len(p.bond_elements)
            print(p.tag, p.type.name, p.type.elem, p.type.desc, p.bond_elements)

    def apply_forcefield(self, f, charges='default', set_box=True, box_padding=10,
                         update_ptypes=False, skip_ptypes=False):
        """pysimm.system.System.apply_forcefield

        Applies force field data to System based on typing rules defined in Forcefield object f

        Args:
            f: pysimm.forcefield.Forcefield object from which new force field type will be retrieved
            charges: type of charges ot be applied default='default'
            set_box: Update simulation box information based on particle positions default=True
            box_padding: Add padding to simulation box if updating dimensions default=10 (Angstroms)
            update_ptypes: If True, update particle types based on current ParticleType names default=False
            skip_ptypes: if True, do not change particle types

        Returns:
            None
        """

        self.ff_class = f.ff_class
        self.ff_name = f.ff_name
        if update_ptypes:
            self.update_particle_types_from_forcefield(f)
            skip_ptypes = True
        if not skip_ptypes:
            f.assign_ptypes(self)
        if self.bonds.count > 0:
            f.assign_btypes(self)
        f.assign_atypes(self)
        f.assign_dtypes(self)
        f.assign_itypes(self)

        if charges:
            f.assign_charges(self, charges=charges)

        if set_box:
            self.set_box(box_padding)

    def apply_charges(self, f, charges='default'):
        f.assign_charges(self, charges=charges)

    def write_amber(self, prmtop='prmtop', inpcrd='inpcrd', **kwargs):
        """pysimm.system.System.write_amber

        *** NOT FINISHED ***

        """
        periodic = kwargs.get('periodic') or 1

        self.set_atomic_numbers()
        self.set_excluded_particles()

        nbonh = 0
        for b in self.bonds:
            if b.a.type.elem == 'H' or b.b.type.elem == 'H':
                nbonh += 1
            if periodic is None:
                if b.distance() > self.dim.dx/2 or b.distance() > self.dim.dy/2 or b.distance() > self.dim.dz/2:
                    periodic = True

        if periodic is not None:
            periodic = 1
        else:
            periodic = 0

        ntheth = 0
        for a in self.angles:
            if a.a.type.elem == 'H' or a.b.type.elem == 'H' or a.c.type.elem == 'H':
                ntheth += 1

        nphih = 0
        for d in self.dihedrals:
            if d.a.type.elem == 'H' or d.b.type.elem == 'H' or d.c.type.elem == 'H' or d.d.type.elem == 'H':
                nphih += 1

        out_prmtop = open(prmtop, 'w+')
        out_inpcrd = open(inpcrd, 'w+')

        out_prmtop.write('%VERSION  WRITTEN BY PYSIMM\n')
        out_prmtop.write('%FLAG TITLE\n')
        out_prmtop.write('%FORMAT(20a4)\n')
        out_prmtop.write('%s\n' % self.name)
        out_prmtop.write('%FLAG POINTERS\n')
        out_prmtop.write('%FORMAT(10I8)\n')
        out_prmtop.write('%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (self.particles.count, self.particle_types.count,
                                                               nbonh, (self.bonds.count - nbonh),
                                                               ntheth, (self.angles.count - ntheth),
                                                               nphih, (self.dihedrals.count - nphih), 0, 0))
        out_prmtop.write('%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (0, self.molecules.count, 0, 0, 0, self.bond_types.count,
                                                               self.angles_types.count, self.dihedral_types.count, 0,
                                                               0))
        out_prmtop.write('%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (0, 0, 0, 0, 0, 0, 0, periodic, 0, 0))
        out_prmtop.write('%8d\n' % 0)

        out_prmtop.write('%FLAG ATOM_NAME\n')
        out_prmtop.write('%FORMAT(20a4)\n')
        for i in range(1, self.particles.count + 1):
            if (i != 1 and i % 20 == 0) or i == self.particles.count:
                out_prmtop.write('%s\n' % self.particles[i].type.name[:4])
            else:
                out_prmtop.write('%s' % self.particles[i].type.name[:4])

        out_prmtop.write('%FLAG CHARGE\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.particles.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.particles.count:
                out_prmtop.write('%.8e\n' % self.particles[i].charge)
            else:
                out_prmtop.write('%.8e' % self.particles[i].charge)

        out_prmtop.write('%FLAG ATOMIC_NUMBER\n')
        out_prmtop.write('%FORMAT(10I8)\n')
        for i in range(1, self.particles.count + 1):
            if (i != 1 and i % 10 == 0) or i == self.particles.count:
                out_prmtop.write('%s\n' % self.particles[i].type.atomic_number)
            else:
                out_prmtop.write('%s' % self.particles[i].type.atomic_number)

        out_prmtop.write('%FLAG MASS\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.particles.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.particles.count:
                out_prmtop.write('%.8e\n' % self.particles[i].type.mass)
            else:
                out_prmtop.write('%.8e' % self.particles[i].type.mass)

        out_prmtop.write('%FLAG ATOM_TYPE_INDEX\n')
        out_prmtop.write('%FORMAT(10I8)\n')
        for i in range(1, self.particles.count + 1):
            if (i != 1 and i % 10 == 0) or i == self.particles.count:
                out_prmtop.write('%s\n' % self.particles[i].type.tag)
            else:
                out_prmtop.write('%s' % self.particles[i].type.tag)

        out_prmtop.write('%FLAG NUMBER_EXCLUDED_ATOMS\n')
        out_prmtop.write('%FORMAT(10I8)\n')
        for i in range(1, self.particles.count + 1):
            excluded = self.particles[i].excluded_particles.count or 1
            if (i != 1 and i % 10 == 0) or i == self.particles.count:
                out_prmtop.write('%s\n' % excluded)
            else:
                out_prmtop.write('%s' % excluded)

        out_prmtop.write('%FLAG NONBONDED_PARM_INDEX\n')
        out_prmtop.write('%FORMAT(10I8)\n')
        count = 0
        for i in range(1, self.particles.count + 1):
            for j in range(1, self.particles.count + 1):
                count += 1
                if i <= j:
                    index = self.particle_types.count * (self.particles[i].type.tag-1) + self.particles[j].type.tag
                else:
                    index = self.particle_types.count * (self.particles[j].type.tag-1) + self.particles[i].type.tag
                if (count != 1 and count % 10 == 0) or count == self.particle_types.count*self.particle_types.count:
                    out_prmtop.write('%s\n' % index)
                else:
                    out_prmtop.write('%s' % index)

        out_prmtop.write('%FLAG NONBONDED_PARM_INDEX\n')
        out_prmtop.write('%FORMAT(10I8)\n')
        count = 0
        for i in range(1, self.particles.count + 1):
            for j in range(1, self.particles.count + 1):
                count += 1
                if i <= j:
                    index = self.particle_types.count * (self.particles[i].type.tag-1) + self.particles[j].type.tag
                else:
                    index = self.particle_types.count * (self.particles[j].type.tag-1) + self.particles[i].type.tag
                if (count != 1 and count % 10 == 0) or count == self.particle_types.count*self.particle_types.count:
                    out_prmtop.write('%s\n' % index)
                else:
                    out_prmtop.write('%s' % index)

        out_prmtop.write('%FLAG RESIDUE_LABEL\n')
        out_prmtop.write('%FORMAT(20a4)\n')
        for i in range(1, self.molecules.count + 1):
            if (i != 1 and i % 20 == 0) or i == self.molecules.count:
                out_prmtop.write('%s\n' % self.molecules[i].name[:4])
            else:
                out_prmtop.write('%s' % self.molecules[i].name[:4])

        out_prmtop.write('%FLAG RESIDUE_POINTER\n')
        out_prmtop.write('%FORMAT(10I8)\n')
        for i in range(1, self.molecules.count + 1):
            if (i != 1 and i % 10 == 0) or i == self.molecules.count:
                for p in self.molecules[i].particles:
                    out_prmtop.write('%s\n' % p.tag)
                    break
            else:
                for p in self.molecules[i].particles:
                    out_prmtop.write('%s' % p.tag)
                    break

        out_prmtop.write('%FLAG BOND_FORCE_CONSTANT\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.bond_types.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.bond_types.count:
                out_prmtop.write('%.8e\n' % (self.bond_types[i].k*2))
            else:
                out_prmtop.write('%.8e' % (self.bond_types[i].k*2))

        out_prmtop.write('%FLAG BOND_EQUIL_VALUE\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.bond_types.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.bond_types.count:
                out_prmtop.write('%.8e\n' % self.bond_types[i].r0)
            else:
                out_prmtop.write('%.8e' % self.bond_types[i].r0)

        out_prmtop.write('%FLAG ANGLE_FORCE_CONSTANT\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.angle_types.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.angle_types.count:
                out_prmtop.write('%.8e\n' % (self.angle_types[i].k*2))
            else:
                out_prmtop.write('%.8e' % (self.angle_types[i].k*2))

        out_prmtop.write('%FLAG ANGLE_EQUIL_VALUE\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.angle_types.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.angle_types.count:
                out_prmtop.write('%.8e\n' % (self.angle_types[i].theta0*pi/180.))
            else:
                out_prmtop.write('%.8e' % (self.angle_types[i].theta0*pi/180.))

        out_prmtop.write('%FLAG DIHEDRAL_FORCE_CONSTANT\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.dihedral_types.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.dihedral_types.count:
                out_prmtop.write('%.8e\n' % self.dihedral_types[i].k)
            else:
                out_prmtop.write('%.8e' % self.dihedral_types[i].k)

        out_prmtop.write('%FLAG DIHEDRAL_PERIODICITY\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.dihedral_types.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.dihedral_types.count:
                out_prmtop.write('%.8e\n' % self.dihedral_types[i].n)
            else:
                out_prmtop.write('%.8e' % self.dihedral_types[i].n)

        out_prmtop.write('%FLAG DIHEDRAL_PHASE\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.dihedral_types.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.dihedral_types.count:
                out_prmtop.write('%.8e\n' % acos(self.dihedral_types[i].d))
            else:
                out_prmtop.write('%.8e' % acos(self.dihedral_types[i].d))

        out_prmtop.write('%FLAG SCEE_SCALE_FACTOR\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.dihedral_types.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.dihedral_types.count:
                out_prmtop.write('%.8e\n' % 1.2)
            else:
                out_prmtop.write('%.8e' % 1.2)

        out_prmtop.write('%FLAG SCNB_SCALE_FACTOR\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        for i in range(1, self.dihedral_types.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.dihedral_types.count:
                out_prmtop.write('%.8e\n' % 2)
            else:
                out_prmtop.write('%.8e' % 2)

        out_prmtop.write('%FLAG LENNARD_JONES_ACOEF\n')
        out_prmtop.write('%FORMAT(5E16.8)\n')
        count = 0
        for i in range(1, self.particle_types.count + 1):
            for j in range(1, self.particle_types.count + 1):
                count += 1
                if i <= j:
                    index = self.particle_types.count * (self.particles[i].type.tag-1) + self.particles[j].type.tag
                else:
                    index = self.particle_types.count * (self.particles[j].type.tag-1) + self.particles[i].type.tag
                if (count != 1 and count % 10 == 0) or count == self.particle_types.count*self.particle_types.count:
                    out_prmtop.write('%s\n' % index)
                else:
                    out_prmtop.write('%s' % index)
        for i in range(1, self.particle_types.count + 1):
            if (i != 1 and i % 5 == 0) or i == self.particle_types.count:
                out_prmtop.write('%.8e\n' % self.particle_types[i].epsilon)
            else:
                out_prmtop.write('%.8e' % self.particle_types[i].epsilon)

    def write_lammps(self, out_data, **kwargs):
        """pysimm.system.System.write_lammps

        Write System data formatted for LAMMPS

        Args:
            out_data: where to write data, file name or 'string'

        Returns:
            None or string of data file if out_data='string'
        """
        empty = kwargs.get('empty')

        if out_data == 'string':
            out_file = StringIO()
        else:
            out_file = open(out_data, 'w+')

        if empty:
            out_file.write('%s\n\n' % self.name)
            out_file.write('%s atoms\n' % 0)
            out_file.write('%s bonds\n' % 0)
            out_file.write('%s angles\n' % 0)
            out_file.write('%s dihedrals\n' % 0)
            out_file.write('%s impropers\n' % 0)
        else:
            out_file.write('%s\n\n' % self.name)
            out_file.write('%s atoms\n' % self.particles.count)
            out_file.write('%s bonds\n' % self.bonds.count)
            out_file.write('%s angles\n' % self.angles.count)
            out_file.write('%s dihedrals\n' % self.dihedrals.count)
            out_file.write('%s impropers\n' % self.impropers.count)

        out_file.write('\n')

        out_file.write('%s atom types\n' % self.particle_types.count)
        if self.bond_types.count > 0:
            out_file.write('%s bond types\n' % self.bond_types.count)
        if self.angle_types.count > 0:
            out_file.write('%s angle types\n' % self.angle_types.count)
        if self.dihedral_types.count > 0:
            out_file.write('%s dihedral types\n' % self.dihedral_types.count)
        if self.improper_types.count > 0:
            out_file.write('%s improper types\n' % self.improper_types.count)

        out_file.write('\n')

        out_file.write('%f %f xlo xhi\n' % (self.dim.xlo, self.dim.xhi))
        out_file.write('%f %f ylo yhi\n' % (self.dim.ylo, self.dim.yhi))
        out_file.write('%f %f zlo zhi\n' % (self.dim.zlo, self.dim.zhi))

        out_file.write('\n')

        if self.particle_types.count > 0:
            out_file.write('Masses\n\n')
            for pt in self.particle_types:
                if not pt.mass:
                    error_print('error: some particle types do not have masses')
                    return
                out_file.write('%4d\t%s\t# %s\n' % (pt.tag, pt.mass, pt.name))
            out_file.write('\n')

        if self.particle_types.count > 0:
            out_file.write('Pair Coeffs\n\n')
            for pt in self.particle_types:
                if (self.pair_style and (self.pair_style.startswith('lj') or
                        self.pair_style.startswith('class2')) and
                        pt.sigma is not None and pt.epsilon is not None):
                    out_file.write('%4d\t%s\t%s\t# %s\n'
                                   % (pt.tag, pt.epsilon, pt.sigma, pt.name))
                elif (self.pair_style and self.pair_style.startswith('buck') and
                        pt.a is not None and pt.rho is not None and pt.c is not None):
                    out_file.write('%4d\t%s\t%s\t%s\t# %s\n'
                                   % (pt.tag, pt.a, pt.rho, pt.c, pt.name))
                elif not self.pair_style and pt.sigma is not None and pt.epsilon is not None:
                    out_file.write('%4d\t%s\t%s\t# %s\n'
                                   % (pt.tag, pt.epsilon, pt.sigma, pt.name))
                elif not self.pair_style and pt.a is not None and pt.rho is not None and pt.c is not None:
                    out_file.write('%4d\t%s\t%s\t%s\t# %s\n'
                                   % (pt.tag, pt.a, pt.rho, pt.c, pt.name))
                else:
                    error_print('error: cannot understand your pair style')
                    return
            out_file.write('\n')

        if self.bond_types.count > 0:
            out_file.write('Bond Coeffs\n\n')
            for b in self.bond_types:
                if self.bond_style == 'harmonic' or self.ff_class == '1':
                    out_file.write('%4d\t%s\t%s\t# %s\n'
                                   % (b.tag, b.k, b.r0, b.name))
                elif self.bond_style == 'class2' or self.ff_class == '2':
                    out_file.write('%4d\t%s\t%s\t%s\t%s\t# %s\n'
                                   % (b.tag, b.r0, b.k2, b.k3, b.k4, b.name))
                else:
                    error_print('error: cannot understand your bond style')
            out_file.write('\n')

        if self.angle_types.count > 0:
            out_file.write('Angle Coeffs\n\n')
            for a in self.angle_types:
                if self.angle_style == 'harmonic' or self.ff_class == '1':
                    out_file.write('%4d\t%s\t'
                                   '%s\t# %s\n'
                                   % (a.tag, a.k,
                                      a.theta0, a.name))
                elif self.angle_style == 'class2' or self.ff_class == '2':
                    out_file.write('%4d\t%s\t'
                                   '%s\t%s\t%s\t# %s\n'
                                   % (a.tag, a.theta0,
                                      a.k2, a.k3, a.k4, a.name))
            out_file.write('\n')

        if (self.angle_types.count > 0 and (self.ff_class == '2' or
                                            self.angle_style == 'class2')):
            out_file.write('BondBond Coeffs\n\n')
            for a in self.angle_types:
                if not a.m:
                    a.m = 0.0
                    if not a.r1:
                        a.r1 = 0.0
                    if not a.r2:
                        a.r2 = 0.0
                out_file.write('%4d\t%s\t%s\t%s\t# %s\n'
                               % (a.tag, a.m, a.r1, a.r2, a.name))
            out_file.write('\n')
            out_file.write('BondAngle Coeffs\n\n')
            for a in self.angle_types:
                if not a.n1:
                    a.n1 = 0.0
                    if not a.r1:
                        a.r1 = 0.0
                    if not a.r2:
                        a.r2 = 0.0
                if not a.n2:
                    a.n2 = 0.0
                out_file.write('%4d\t%s\t%s\t%s\t%s\t# %s\n'
                               % (a.tag, a.n1, a.n2, a.r1, a.r2, a.name))
            out_file.write('\n')

        if self.dihedral_types.count > 0:
            out_file.write('Dihedral Coeffs\n\n')
            for dt in self.dihedral_types:
                if self.dihedral_style == 'harmonic' or self.ff_class == '1':
                    out_file.write('%4d\t%s\t%2s\t%s\t# %s\n'
                                   % (dt.tag, dt.k, dt.d, dt.n, dt.name))
                elif self.dihedral_style == 'fourier':
                    dt_str = '{}'.format(dt.m)
                    for k, d, n in zip(dt.k, dt.d. dt.n):
                        dt_str += '\t{}\t{}\t{}'.format(k, d, n)
                    dt_str += '\t# {}\n'.format(dt.name)
                    out_file.write(dt_str)
                elif self.dihedral_style == 'class2' or self.ff_class == '2':
                    out_file.write('%4d\t%s\t%s\t%s\t%s\t%s\t%s\t# %s\n'
                                   % (dt.tag, dt.k1, dt.phi1, dt.k2, dt.phi2, dt.k3,
                                      dt.phi3, dt.name))
            out_file.write('\n')

        if self.dihedral_types.count > 0 and (self.ff_class == '2' or
                                        self.dihedral_style == 'class2'):
            out_file.write('MiddleBondTorsion Coeffs\n\n')
            for d in self.dihedral_types:
                if not d.a1:
                    d.a1 = 0.0
                    if not d.r2:
                        d.r2 = 0.0
                if not d.a2:
                    d.a2 = 0.0
                    if not d.r2:
                        d.r2 = 0.0
                if not d.a3:
                    d.a3 = 0.0
                    if not d.r2:
                        d.r2 = 0.0
                out_file.write('%4d\t%s\t%s\t%s\t%s\t# %s\n'
                               % (d.tag, d.a1, d.a2, d.a3, d.r2, d.name))
            out_file.write('\n')
            out_file.write('EndBondTorsion Coeffs\n\n')
            for d in self.dihedral_types:
                if not d.b1:
                    d.b1 = 0.0
                    if not d.r1:
                        d.r1 = 0.0
                    if not d.r3:
                        d.r3 = 0.0
                if not d.b2:
                    d.b2 = 0.0
                    if not d.r1:
                        d.r1 = 0.0
                    if not d.r3:
                        d.r3 = 0.0
                if not d.b3:
                    d.b3 = 0.0
                    if not d.r1:
                        d.r1 = 0.0
                    if not d.r3:
                        d.r3 = 0.0
                if not d.c1:
                    d.c1 = 0.0
                    if not d.r1:
                        d.r1 = 0.0
                    if not d.r3:
                        d.r3 = 0.0
                if not d.c2:
                    d.c2 = 0.0
                    if not d.r1:
                        d.r1 = 0.0
                    if not d.r3:
                        d.r3 = 0.0
                if not d.c3:
                    d.c3 = 0.0
                    if not d.r1:
                        d.r1 = 0.0
                    if not d.r3:
                        d.r3 = 0.0
                out_file.write('%4d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t# %s\n'
                               % (d.tag,
                                  d.b1, d.b2, d.b3,
                                  d.c1, d.c2, d.c3,
                                  d.r1, d.r3,
                                  d.name))
            out_file.write('\n')
            out_file.write('AngleTorsion Coeffs\n\n')
            for d in self.dihedral_types:
                if not d.d1:
                    d.d1 = 0.0
                    if not d.theta1:
                        d.theta1 = 0.0
                    if not d.theta2:
                        d.theta2 = 0.0
                if not d.d2:
                    d.d2 = 0.0
                    if not d.theta1:
                        d.theta1 = 0.0
                    if not d.theta2:
                        d.theta2 = 0.0
                if not d.d3:
                    d.d3 = 0.0
                    if not d.theta1:
                        d.theta1 = 0.0
                    if not d.theta2:
                        d.theta2 = 0.0
                if not d.e1:
                    d.e1 = 0.0
                    if not d.theta1:
                        d.theta1 = 0.0
                    if not d.theta2:
                        d.theta2 = 0.0
                if not d.e2:
                    d.e2 = 0.0
                    if not d.theta1:
                        d.theta1 = 0.0
                    if not d.theta2:
                        d.theta2 = 0.0
                if not d.e3:
                    d.e3 = 0.0
                    if not d.theta1:
                        d.theta1 = 0.0
                    if not d.theta2:
                        d.theta2 = 0.0
                out_file.write('%4d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t# %s\n'
                               % (d.tag,
                                  d.d1, d.d2, d.d3,
                                  d.e1, d.e2, d.e3,
                                  d.theta1, d.theta2,
                                  d.name))
            out_file.write('\n')
            out_file.write('AngleAngleTorsion Coeffs\n\n')
            for d in self.dihedral_types:
                if not d.m:
                    d.m = 0.0
                    if not d.theta1:
                        d.theta1 = 0.0
                    if not d.theta2:
                        d.theta2 = 0.0
                out_file.write('%4d\t%s\t%s\t%s\t# %s\n'
                               % (d.tag,
                                  d.m,
                                  d.theta1, d.theta2,
                                  d.name))
            out_file.write('\n')
            out_file.write('BondBond13 Coeffs\n\n')
            for d in self.dihedral_types:
                if not d.n_class2:
                    d.n_class2 = 0.0
                    if not d.r1:
                        d.r1 = 0.0
                    if not d.r3:
                        d.r3 = 0.0
                out_file.write('%4d\t%s\t%s\t%s\t# %s\n'
                               % (d.tag,
                                  d.n_class2,
                                  d.r1, d.r3,
                                  d.name))
            out_file.write('\n')

        if self.improper_types.count > 0:
            out_file.write('Improper Coeffs\n\n')
            for i in self.improper_types:
                if self.improper_style = 'harmonic' or self.improper_style =='class2':
                    if not i.k:
                        i.k = 0.0
                    if not i.x0:
                        i.x0 = 0.0
                    out_file.write('%4d\t%s\t%s\t# %s\n'
                                   % (i.tag, i.k, i.x0, i.name))
                elif self.improper_style = 'cvff':
                    out_file.write('%4d\t%s\t%s\t%s\t# %s\n'
                                   % (i.tag, i.k, i.d, i.n, i.name))
            out_file.write('\n')

        if self.improper_types.count > 0 and (self.ff_class == '2' or
                                              self.improper_style == 'class2'):
            out_file.write('AngleAngle Coeffs\n\n')
            for i in self.improper_types:
                if not i.m1:
                    i.m1 = 0.0
                    if not i.theta1:
                        i.theta1 = 0.0
                    if not i.theta2:
                        i.theta2 = 0.0
                    if not i.theta3:
                        i.theta3 = 0.0
                if not i.m2:
                    i.m2 = 0.0
                    if not i.theta1:
                        i.theta1 = 0.0
                    if not i.theta2:
                        i.theta2 = 0.0
                    if not i.theta3:
                        i.theta3 = 0.0
                if not i.m3:
                    i.m3 = 0.0
                    if not i.theta1:
                        i.theta1 = 0.0
                    if not i.theta2:
                        i.theta2 = 0.0
                    if not i.theta3:
                        i.theta3 = 0.0
                out_file.write('%4d\t%s\t%s\t%s\t%s\t%s\t%s\t# %s\n'
                               % (i.tag,
                                  i.m1, i.m2, i.m3,
                                  i.theta1, i.theta2, i.theta3,
                                  i.name))
            out_file.write('\n')

        if self.particles.count > 0 and not empty:
            out_file.write('Atoms\n\n')
            for p in self.particles:
                if not p.molecule:
                    p.molecule = Item()
                    p.molecule.tag = 1
                if not p.charge:
                    p.charge = 0
                if isinstance(p.molecule, int):
                    out_file.write('%4d\t%d\t%d\t%s\t%s\t%s\t%s\n'
                                   % (p.tag, p.molecule, p.type.tag, p.charge,
                                      p.x, p.y, p.z))
                else:
                    out_file.write('%4d\t%d\t%d\t%s\t%s\t%s\t%s\n'
                                   % (p.tag, p.molecule.tag, p.type.tag, p.charge,
                                      p.x, p.y, p.z))
            out_file.write('\n')

            out_file.write('Velocities\n\n')
            for p in self.particles:
                if not p.vx:
                    p.vx = 0.
                if not p.vy:
                    p.vy = 0.
                if not p.vz:
                    p.vz = 0
                out_file.write('%4d\t%s\t%s\t%s\n' % (p.tag, p.vx, p.vy, p.vz))
            out_file.write('\n')

        if self.bonds.count > 0 and not empty:
            out_file.write('Bonds\n\n')
            for b in self.bonds:
                out_file.write('%4d\t%d\t%d\t%d\n'
                               % (b.tag, b.type.tag, b.a.tag, b.b.tag))
            out_file.write('\n')

        if self.angles.count > 0 and not empty:
            out_file.write('Angles\n\n')
            for a in self.angles:
                out_file.write('%4d\t%d\t%d\t%d\t%d\n'
                               % (a.tag, a.type.tag, a.a.tag, a.b.tag, a.c.tag))
            out_file.write('\n')

        if self.dihedrals.count > 0 and not empty:
            out_file.write('Dihedrals\n\n')
            for d in self.dihedrals:
                out_file.write('%4d\t%d\t%d\t%d\t%d\t%d\n'
                               % (d.tag, d.type.tag,
                                  d.a.tag, d.b.tag, d.c.tag, d.d.tag))
            out_file.write('\n')

        if self.impropers.count > 0 and not empty:
            out_file.write('Impropers\n\n')
            for i in self.impropers:
                if self.ff_class == '2' or self.improper_style == 'class2':
                    out_file.write('%4d\t%d\t%d\t%d\t%d\t%d\n'
                                   % (i.tag, i.type.tag,
                                      i.b.tag, i.a.tag, i.c.tag, i.d.tag))
                else:
                    out_file.write('%4d\t%d\t%d\t%d\t%d\t%d\n'
                                   % (i.tag, i.type.tag,
                                      i.a.tag, i.b.tag, i.c.tag, i.d.tag))
            out_file.write('\n')

        if out_data == 'string':
            s = out_file.getvalue()
            out_file.close()
            return s
        else:
            out_file.close()

    def write_hoomd(self, outfile='data.xml'):
        """pysimm.system.System.write_hoomd

        Write System data formatted for hoomd

        Args:
            outfile: file name to write data

        Returns:
            None
        """
        hoomd = Et.Element('hoomd_xml')
        tree = Et.ElementTree(hoomd)

        config = Et.SubElement(hoomd, 'configuration')

        box = Et.SubElement(config, 'box')
        box.set('lx', str((self.dim.xhi * 0.1) - (self.dim.xlo * 0.1)))
        box.set('ly', str((self.dim.yhi * 0.1) - (self.dim.ylo * 0.1)))
        box.set('lz', str((self.dim.zhi * 0.1) - (self.dim.zlo * 0.1)))

        position = Et.SubElement(config, 'position')
        position.set('num', str(self.particles.count))
        position.text = ''
        for p in self.particles:
            position.text += ('%s %s %s\n' % (p.x * 0.1, p.y * 0.1, p.z * 0.1))

        mass = Et.SubElement(config, 'mass')
        mass.set('num', str(self.particles.count))
        mass.text = ''
        for p in self.particles:
            mass.text += '%s\n' % p.type.mass

        type_ = Et.SubElement(config, 'type')
        type_.set('num', str(self.particles.count))
        type_.text = ''
        for p in self.particles:
            type_.text += '%s\n' % p.type.name

        if self.bonds.count > 0:
            bond = Et.SubElement(config, 'bond')
            bond.set('num', str(self.bonds.count))
            bond.text = ''
            for b in self.bonds:
                bond.text += '%s %s %s\n' % (b.type.name, b.a.tag - 1,
                                             b.b.tag - 1)

        if self.angles.count > 0:
            angle = Et.SubElement(config, 'angle')
            angle.set('num', str(self.angles.count))
            angle.text = ''
            for a in self.angles:
                angle.text += '%s %s %s %s\n' % (a.type.name, a.a.tag - 1,
                                                 a.b.tag - 1, a.c.tag - 1)

        if self.dihedrals.count > 0:
            dihedral = Et.SubElement(config, 'dihedral')
            dihedral.set('num', str(self.dihedrals.count))
            dihedral.text = ''
            for d in self.dihedrals:
                dihedral.text += '%s %s %s %s %s\n' % (d.type.name,
                                                       d.a.tag - 1,
                                                       d.b.tag - 1,
                                                       d.c.tag - 1,
                                                       d.d.tag - 1)

        if self.impropers.count > 0:
            improper = Et.SubElement(config, 'improper')
            improper.set('num', str(self.impropers.count))
            improper.text = ''
            for i in self.impropers:
                improper.text += '%s %s %s %s %s\n' % (i.type.name,
                                                       i.a.tag - 1,
                                                       i.b.tag - 1,
                                                       i.c.tag - 1,
                                                       i.d.tag - 1)

        tree.write(outfile)

    def write_xyz(self, outfile='data.xyz', **kwargs):
        """pysimm.system.System.write_xyz

        Write System data in xyz format

        Args:
            outfile: where to write data, file name or 'string'

        Returns:
            None or string of data file if out_data='string'
        """
        elem = kwargs.get('elem') if kwargs.get('elem') is not None else True
        append = kwargs.get('append')
        if outfile == 'string':
            out = StringIO()
        else:
            if append:
                out = open(outfile, 'a')
            else:
                out = open(outfile, 'w')

        out.write('%s\n' % self.particles.count)
        out.write('xyz file written from pySIMM system module\n')
        for p in self.particles:
            if elem and p.type and p.type.elem is not None:
                out.write('%s %s %s %s\n' % (p.type.elem, p.x, p.y, p.z))
            elif elem and p.elem is not None:
                out.write('%s %s %s %s\n' % (p.elem, p.x, p.y, p.z))
            else:
                out.write('%s %s %s %s\n' % (p.type.tag, p.x, p.y, p.z))
        if outfile == 'string':
            s = out.getvalue()
            out.close()
            return s
        else:
            out.close()

    def write_chemdoodle_json(self, outfile, **kwargs):
        """pysimm.system.System.write_chemdoodle_json

        Write System data in chemdoodle json format

        Args:
            outfile: where to write data, file name or 'string'

        Returns:
            None or string of data file if out_data='string'
        """

        atoms = []
        bonds = []

        for p in self.particles:
            if p.type and p.type.elem:
                atoms.append({"x": p.x, "y": p.y, "z": p.z, "l": p.type.elem, "i": p.type.name, "c": p.charge})
            elif p.elem and p.type:
                atoms.append({"x": p.x, "y": p.y, "z": p.z, "l": p.elem, "i": p.type.name, "c": p.charge})
            elif p.elem:
                atoms.append({"x": p.x, "y": p.y, "z": p.z, "l": p.elem})
            else:
                atoms.append({"x": p.x, "y": p.y, "z": p.z, "i": p.type.name, "c": p.charge})

        for b in self.bonds:
            if b.order:
                bonds.append({"b": b.a.tag-1, "e": b.b.tag-1, "o": b.order})
            else:
                bonds.append({"b": b.a.tag-1, "e": b.b.tag-1})

        j = {"a": atoms, "b": bonds}

        if outfile == 'string':
            out = StringIO()
        else:
            out = open(outfile, 'w+')

        out.write(json.dumps(j))

        if outfile == 'string':
            s = out.getvalue()
            out.close()
            return s
        else:
            out.close()

    def write_mol(self, outfile='data.mol'):
        """pysimm.system.System.write_mol

        Write System data in mol format

        Args:
            outfile: where to write data, file name or 'string'

        Returns:
            None or string of data file if out_data='string'
        """
        if outfile == 'string':
            out = StringIO()
        else:
            out = open(outfile, 'w+')

        out.write('system\n')
        out.write('written using pySIMM system module\n\n')
        out.write('%s\t%s\n' % (self.particles.count, self.bonds.count))
        for p in self.particles:
            if not p.charge:
                p.charge = 0.0
            if p.type and p.type.elem:
                out.write('%10.4f%10.4f%10.4f %s 0 %10.4f\n'
                          % (p.x, p.y, p.z, '{0: >3}'.format(p.type.elem),
                             p.charge))
            elif p.elem:
                out.write('%10.4f%10.4f%10.4f %s 0 %10.4f\n'
                          % (p.x, p.y, p.z, '{0: >3}'.format(p.elem),
                             p.charge))
            elif p.type:
                out.write('%10.4f%10.4f%10.4f %s 0 %10.4f\n'
                          % (p.x, p.y, p.z, '{0: >3}'.format(p.type.tag),
                             p.charge))

        for b in self.bonds:
            if b.order:
                out.write('%s\t%s\t%s\t%s\t%s\t%s\n'
                          % (b.a.tag, b.b.tag, b.order, 0, 0, 0))
            else:
                out.write('%s\t%s\t%s\t%s\t%s\t%s\n'
                          % (b.a.tag, b.b.tag, 1, 0, 0, 0))
        out.write('M END')
        if outfile == 'string':
            s = out.getvalue()
            out.close()
            return s
        else:
            out.close()

    def write_pdb(self, outfile='data.pdb'):
        """pysimm.system.System.write_pdb

        Write System data in pdb format

        Args:
            outfile: where to write data, file name or 'string'

        Returns:
            None or string of data file if out_data='string'
        """
        if outfile == 'string':
            out = StringIO()
        else:
            out = open(outfile, 'w+')

        out.write('{0: <10}pdb written using pySIMM system module\n'
                  .format('HEADER'))
        for p in self.particles:
            out.write('{0: <6}{1: >5} {2: >4} RES {3}     '
                      '{4: >8.3f}{5: >8.3f}{6: >8.3f}{7: >24}{8: >2}\n'
                      .format('ATOM', p.tag, p.type.name, p.molecule.tag,
                              p.x, p.y, p.z, '', p.type.elem))
        for p in self.particles:
            if p.bonds:
                out.write('{0: <6}{1: >5}'
                          .format('CONECT', p.tag))
                for t in sorted([x.a.tag if p is x.b else x.b.tag for x in
                                 p.bonds]):
                    out.write('{0: >5}'.format(t))
                out.write('\n')

        if outfile == 'string':
            s = out.getvalue()
            out.close()
            return s
        else:
            out.close()

    def write_yaml(self, file_):
        """pysimm.system.System.write_yaml

        Write System data in yaml format

        Args:
            outfile: file name to write data

        Returns:
            None
        """
        n = self.copy()

        s = vars(n)
        for k, v in s.items():
            if isinstance(v, ItemContainer):
                s[k] = vars(v)
                for k_, v_ in s[k].items():
                    if k_ == '_dict':
                        for t, i in v_.items():
                            s[k][k_][t] = vars(i)
                            for key, value in s[k][k_][t].items():
                                if isinstance(value, ItemContainer) or (isinstance(value, list) and
                                                                        value and isinstance(value[0], Item)):
                                    s[k][k_][t][key] = [x.tag for x in value]
                                elif isinstance(value, Item) or isinstance(value, System) and value.tag:
                                    s[k][k_][t][key] = value.tag
            elif isinstance(v, Item):
                s[k] = vars(v)

        if file_ == 'string':
            f = StringIO()
            f.write(json.dumps(s, indent=4, separators=(',', ': ')))
            yaml_ = f.getvalue()
            f.close()
            return yaml_

        with file(file_, 'w') as f:
            f.write(json.dumps(s, indent=4, separators=(',', ': ')))

    def consolidate_types(self):
        """pysimm.system.System.consolidate_types

        Removes duplicate types and reassigns references

        Args:
            None

        Returns:
            None
        """
        for pt in self.particle_types:
            for dup in self.particle_types:
                if pt is not dup and pt.name == dup.name:
                    for p in self.particles:
                        if p.type == dup:
                            p.type = pt
                    self.particle_types.remove(dup.tag)

        for bt in self.bond_types:
            for dup in self.bond_types:
                if bt is not dup and bt.name == dup.name:
                    for b in self.bonds:
                        if b.type == dup:
                            b.type = bt
                    self.bond_types.remove(dup.tag)

        for at in self.angle_types:
            for dup in self.angle_types:
                if at is not dup and at.name == dup.name:
                    for a in self.angles:
                        if a.type == dup:
                            a.type = at
                    self.angle_types.remove(dup.tag)

        for dt in self.dihedral_types:
            for dup in self.dihedral_types:
                if dt is not dup and dt.name == dup.name:
                    for d in self.dihedrals:
                        if d.type == dup:
                            d.type = dt
                    self.dihedral_types.remove(dup.tag)

        for it in self.improper_types:
            for dup in self.improper_types:
                if it is not dup and it.name == dup.name:
                    for i in self.impropers:
                        if i.type == dup:
                            i.type = it
                    self.improper_types.remove(dup.tag)

    def set_cog(self):
        """pysimm.system.System.set_cog

        Calculate center of gravity of System and assign to System.cog

        Args:
            None

        Returns:
            None
        """
        self.cog = [0, 0, 0]
        for p in self.particles:
            self.cog[0] += p.x
            self.cog[1] += p.y
            self.cog[2] += p.z
        if self.particles.count:
            self.cog = [c / self.particles.count for c in self.cog]

    def center_at_origin(self):
        """pysimm.system.System.center_at_origin

        Moves particles in system such that new center of gravity is 0, 0, 0

        Args:
            None

        Returns:
            None
        """
        self.set_cog()
        for p in self.particles:
            p.x -= self.cog[0]
            p.y -= self.cog[1]
            p.z -= self.cog[2]
        self.set_cog()

    def set_neighbors(self, cutoff=12.0):
        """pysimm.system.System.set_neighbors

        *** BUGGY - DO NOT USE ***

        """
        ncells = [int(ceil(self.dim.dx/cutoff)), int(ceil(self.dim.dy/cutoff)), int(ceil(self.dim.dz/cutoff))]
        cell_dx = self.dim.dx/ncells[0]
        cell_dy = self.dim.dy/ncells[1]
        cell_dz = self.dim.dz/ncells[2]
        if min(ncells) < 2:
            warning_print('dimension along one axis is less than cutoff - this is not good')
        self.neighbors = [[[ItemContainer() for c_z in range(ncells[2])]
                           for c_y in range(ncells[1])]
                          for c_x in range(ncells[0])]

        for p in self.particles:
            p.cell_id = [int(floor((p.x-self.dim.xlo)/cell_dx)),
                         int(floor((p.y-self.dim.ylo)/cell_dy)),
                         int(floor((p.z-self.dim.zlo)/cell_dz))]
            self.neighbors[p.cell_id[0]][p.cell_id[1]][p.cell_id[2]].add(p)

        for p in self.particles:
            p.neighbors = ItemContainer()
            p_xcell = p.cell_id[0]
            p_ycell = p.cell_id[1]
            p_zcell = p.cell_id[2]

            for x_index in range(p_xcell-1, p_xcell+2):
                for y_index in range(p_ycell-1, p_ycell+2):
                    for z_index in range(p_zcell-1, p_zcell+2):
                        if x_index >= len(self.neighbors):
                            x_index -= len(self.neighbors)
                        if y_index >= len(self.neighbors):
                            y_index -= len(self.neighbors)
                        if z_index >= len(self.neighbors):
                            z_index -= len(self.neighbors)

                        for p_ in self.neighbors[x_index][y_index][z_index]:
                            if p is not p_:
                                p.neighbors.add(p_)

        self.neighbors_check = True

    def guess_bonds(self, heavy_d=2.0, hydrogen_d=1.5, neighbor_cutoff=2.5):
        """pysimm.system.System.guess_bonds

        *** USES BUGGY FUNCTION SET_NEIGHBORS - DO NOT USE ***

        """
        self.bonds.remove('all')
        for p in self.particles:
            p.bonds = ItemContainer()
            p.bonded_to = ItemContainer()

        if not self.neighbors_check:
            self.set_neighbors(cutoff=neighbor_cutoff)

        for p in self.particles:
            for p_ in p.neighbors:
                if p_ in p.bonded_to:
                    continue
                if p.type.elem == 'H' or p_.type.elem == 'H':
                    if calc.pbc_distance(self, p, p_) <= hydrogen_d:
                        new_bond = Bond(a=p, b=p_)
                        self.bonds.add(new_bond)
                        p.bonds.add(new_bond)
                        p_.bonds.add(new_bond)
                        p.bonded_to.add(p_)
                        p_.bonded_to.add(p)
                else:
                    if calc.pbc_distance(self, p, p_) <= heavy_d:
                        new_bond = Bond(a=p, b=p_)
                        self.bonds.add(new_bond)
                        p.bonds.add(new_bond)
                        p_.bonds.add(new_bond)
                        p.bonded_to.add(p_)
                        p_.bonded_to.add(p)

    def guess_hybridization(self, unwrap=True):
        """pysimm.system.System.guess_hybridization

        Simple guessing of hybridization based on element and number of bonds. Cannot detect aromaticity.

        Args:
            unwrap: if True, unwrap System first default=True

        Returns:
            None
        """
        if unwrap:
            self.unwrap()
        self.add_particle_bonding()
        for p in self.particles:
            if not p.type or not p.type.elem:
                error_print('cannot determine hybridization for particle %s because type/element not defined' % p.tag)
                return
            if p.type.elem == 'C':
                if p.bonds.count == 4:
                    p.hyb = 'sp3'
                elif p.bonds.count == 3:
                    p.hyb = 'sp2'
                elif p.bonds.count == 2:
                    p.hyb = 'sp1'
                else:
                    error_print('unexpected number of bonds for particle %s' % p.tag)
            elif p.type.elem == 'O':
                if p.bonds.count == 2:
                    p.hyb = 'sp3'
                elif p.bonds.count == 1:
                    p.hyb = 'sp2'
                else:
                    error_print('unexpected number of bonds for particle %s' % p.tag)
            elif p.type.elem == 'N':
                if p.bonds.count == 1:
                    p.hyb = 'sp1'
                elif p.bonds.count == 2:
                    p.hyb = 'sp2'
                elif p.bonds.count == 3:
                    p.hyb = 'sp3'
                else:
                    error_print('unexpected number of bonds for particle %s' % p.tag)

    def guess_bond_order(self):
        """pysimm.system.System.guess_bond_order

        Simple guessing of bond order based on hybridizations

        Args:
            None

        Returns:
            None
        """
        for b in self.bonds:
            if b.a.type.elem == 'H' or b.b.type.elem == 'H':
                if not b.order:
                    b.order = 1
                else:
                    print(b.order, 1)
            elif b.a.hyb == 'sp2' and b.b.hyb == 'sp2':
                if not b.order:
                    b.order = 2
                else:
                    print(b.order, 2)
            elif b.a.hyb == 'sp1' and b.b.hyb == 'sp1':
                if not b.order:
                    b.order = 3
                else:
                    print(b.order, 3)
            else:
                if not b.order:
                    b.order = 1
                else:
                    print(b.order, 1)

    def set_mass(self):
        """pysimm.system.System.set_mass

        Set total mass of particles in System

        Args:
            None

        Returns:
            None
        """
        self.mass = 0
        for p in self.particles:
            if p.type.mass is None:
                self.mass = 0
                warning_print('Some particles do not have a mass')
                break
            self.mass += p.type.mass

    def set_volume(self):
        """pysimm.system.System.set_volume

        Set volume of System based on Dimension

        Args:
            None

        Returns:
            None
        """
        if self.dim.check():
            self.volume = ((self.dim.xhi - self.dim.xlo) *
                           (self.dim.yhi - self.dim.ylo) *
                           (self.dim.zhi - self.dim.zlo))

    def set_density(self):
        """pysimm.system.System.set_density

        Calculate density of System from mass and volume

        Args:
            None

        Returns:
            None
        """
        self.set_mass()
        self.set_volume()
        if self.mass and self.volume:
            self.density = self.mass / 6.02e23 / self.volume * 1e24

    def set_velocity(self):
        """pysimm.system.System.set_velocity

        Calculate total velocity of particles in System

        Args:
            None

        Returns:
            None
        """
        self.vx = 0.0
        self.vy = 0.0
        self.vz = 0.0
        for p in self.particles:
            if p.vx is None:
                p.vx = 0
            self.vx += p.vx
            if p.vy is None:
                p.vy = 0
            self.vy += p.vy
            if p.vz is None:
                p.vz = 0
            self.vz += p.vz

    def zero_velocity(self):
        """pysimm.system.System.zero_velocity

        Enforce zero shift velocity in system

        Args:
            None

        Returns:
            None
        """
        self.set_velocity()
        shift_x = shift_y = shift_z = 0.0
        if self.vx != 0:
            shift_x = self.vx / self.particles.count
        if self.vy != 0:
            shift_y = self.vy / self.particles.count
        if self.vz != 0:
            shift_z = self.vz / self.particles.count
        if shift_x != 0 or shift_y != 0 or shift_z != 0:
            for p in self.particles:
                p.vx -= shift_x
                p.vy -= shift_y
                p.vz -= shift_z
        self.set_velocity()

    def set_box(self, padding=0., center=True):
        """pysimm.system.System.set_box

        Update System.Dimension with user defined padding

        Args:
            padding: add padding to all sides of box (Angstrom)
            center: if True, place center of box at origin default=True

        Returns:
            None
        """
        if center:
            self.center_at_origin()
        xmin = ymin = zmin = sys.float_info.max
        xmax = ymax = zmax = sys.float_info.min
        for p in self.particles:
            if p.x < xmin:
                xmin = p.x
            if p.x > xmax:
                xmax = p.x
            if p.y < ymin:
                ymin = p.y
            if p.y > ymax:
                ymax = p.y
            if p.z < zmin:
                zmin = p.z
            if p.z > zmax:
                zmax = p.z
        self.dim.xlo = xmin - padding
        self.dim.xhi = xmax + padding
        self.dim.ylo = ymin - padding
        self.dim.yhi = ymax + padding
        self.dim.zlo = zmin - padding
        self.dim.zhi = zmax + padding

        self.dim.dx = self.dim.xhi - self.dim.xlo
        self.dim.dy = self.dim.yhi - self.dim.ylo
        self.dim.dz = self.dim.zhi - self.dim.zlo

    def set_mm_dist(self, molecules=None):
        """pysimm.system.System.set_mm_dist

        Calculate molecular mass distribution (mainly for polymer systems).
        Sets System.mw, System.mn, and System.disperisty

        Args:
            molecules: ItemContainer of molecules to calculate distributions defaul='all'

        Returns:
            None
        """
        if molecules is None or molecules == 'all':
            molecules = self.molecules
        for m in molecules:
            m.set_mass()

        self.mn = 0
        self.mw = 0

        for m in molecules:
            self.mn += m.mass
            self.mw += pow(m.mass, 2)

        self.mw /= self.mn
        self.mn /= molecules.count
        self.dispersity = self.mw / self.mn
        self.pdi = self.mw / self.mn

    def set_frac_free_volume(self, v_void=None):
        """pysimm.system.System.set_frac_free_volume

        Calculates fractional free volume from void volume and bulk density

        Args:
            v_void: void volume if not defined in System.void_volume default=None

        Returns:
            None
        """
        if not v_void and not self.void_volume:
            error_print('Void volume not provided, cannot calculate fractional free volume')
            return
        elif not v_void:
            self.set_density()
            self.frac_free_volume = calc.frac_free_volume(1/self.density, self.void_volume)
        elif not self.void_volume:
            self.set_density()
            self.frac_free_volume = calc.frac_free_volume(1/self.density, v_void)
        if not self.frac_free_volume or self.frac_free_volume < 0:
            self.frac_free_volume = 0.0

    def visualize(self, vis_exec='vmd', **kwargs):
        """pysimm.system.System.visualize

        Visualize system in third party software with given executable. Software must accept pdb or xyz as first
        command line argument.

        Args:
            vis_exec: executable to launch visualization software default='vmd'
            unwrap (optional): if True, unwrap System first default=None
            format (optional): set format default='xyz'

        Returns:
            None
        """

        if not call:
            raise PysimmError('pysimm.system.System.visualize function requires subprocess.call')

        unwrap = kwargs.get('unwrap')
        format = kwargs.get('format') or 'xyz'

        verbose_print(self.dim.dx, self.dim.xlo, self.dim.xhi)
        verbose_print(self.dim.dy, self.dim.ylo, self.dim.yhi)
        verbose_print(self.dim.dz, self.dim.zlo, self.dim.zhi)

        if unwrap:
            self.unwrap()
        if format == 'xyz':
            name_ = 'pysimm_temp.xyz'
            self.write_xyz(name_)
        elif format == 'pdb':
            name_ = 'pysimm_temp.pdb'
            self.write_pdb(name_)

        call('%s %s' % (vis_exec, name_), shell=True)

        os.remove(name_)

    def viz(self, **kwargs):
        self.visualize(vis_exec='vmd', unwrap=False, format='xyz')


class Molecule(System):
    """pysimm.system.Molecule

    Very similar to System, but requires less information
    """
    def __init__(self, **kwargs):
        System.__init__(self, **kwargs)


def join(system1, system2, remove_overlaps=True, max_buffer=6):
    """pysimm.system.join

    *** USES BUGGY FUNCTION SET_NEIGHBORS - DO NOT USE ***

    """
    s1 = system1.copy()
    s2 = system2.copy()
    new_system = replicate([s1, s2], [1, 1], density=None, rand=False, print_insertions=False)
    new_system.dim = s1.dim.copy()
    if not remove_overlaps:
        return new_system
    for m in new_system.molecules[s1.molecules.count+1:]:
        overlap = False
        for mp in m.particles:
            for p in new_system.particles[1:s1.particles.count+1]:
                if calc.pbc_distance(new_system, mp, p) < (mp.type.sigma + p.type.sigma)/2:
                    overlap = True
                    break
            if overlap:
                break
        if overlap:
            new_system.molecules.remove(m.tag, update=False)
            for mp in m.particles:
                new_system.particles.remove(mp.tag, update=False)
    new_system.remove_spare_bonding()
    return new_system


def read_yaml(file_, **kwargs):
    """pysimm.system.read_yaml

    Interprets yaml file and creates pysimm.system.System object

    Args:
        file_: yaml file name

    Returns:
        pysimm.system.System object
    """

    if os.path.isfile(file_):
        dict_ = json.loads(file(file_).read())
    else:
        dict_ = json.loads(file_)

    s = System()

    for k, v in dict_.items():
        if not isinstance(v, dict):
            setattr(s, k, v)

    if isinstance(dict_.get('dim'), dict):
        s.dim = Dimension(**dict_.get('dim'))

    if isinstance(dict_.get('particle_types'), dict):
        s.particle_types = ItemContainer()
        for pt in dict_.get('particle_types').get('_dict').values():
            s.particle_types.add(ParticleType(**pt))

    if isinstance(dict_.get('bond_types'), dict):
        s.bond_types = ItemContainer()
        for bt in dict_.get('bond_types').get('_dict').values():
            s.bond_types.add(BondType(**bt))

    if isinstance(dict_.get('angle_types'), dict):
        s.angle_types = ItemContainer()
        for at in dict_.get('angle_types').get('_dict').values():
            s.angle_types.add(AngleType(**at))

    if isinstance(dict_.get('dihedral_types'), dict):
        s.dihedral_types = ItemContainer()
        for dt in dict_.get('dihedral_types').get('_dict').values():
            s.dihedral_types.add(DihedralType(**dt))

    if isinstance(dict_.get('improper_types'), dict):
        s.improper_types = ItemContainer()
        for it in dict_.get('improper_types').get('_dict').values():
            s.improper_types.add(ImproperType(**it))

    if isinstance(dict_.get('particles'), dict):
        s.particles = ItemContainer()
        for p in dict_.get('particles').get('_dict').values():
            s.particles.add(Particle(**p))

    if isinstance(dict_.get('bonds'), dict):
        s.bonds = ItemContainer()
        for b in dict_.get('bonds').get('_dict').values():
            s.bonds.add(Bond(**b))

    if isinstance(dict_.get('angles'), dict):
        s.angles = ItemContainer()
        for a in dict_.get('angles').get('_dict').values():
            s.angles.add(Angle(**a))

    if isinstance(dict_.get('dihedrals'), dict):
        s.dihedrals = ItemContainer()
        for d in dict_.get('dihedrals').get('_dict').values():
            s.dihedrals.add(Dihedral(**d))

    if isinstance(dict_.get('impropers'), dict):
        s.impropers = ItemContainer()
        for i in dict_.get('impropers').get('_dict').values():
            s.impropers.add(Improper(**i))

    if isinstance(dict_.get('molecules'), dict):
        s.molecules = ItemContainer()
        for m in dict_.get('molecules').get('_dict').values():
            mol = Molecule()
            for k, v in m.items():
                if isinstance(v, list) and not v:
                    setattr(mol, k, ItemContainer())
                else:
                    setattr(mol, k, v)

            particles = [x for x in mol.particles]
            mol.particles = ItemContainer()
            for n in particles:
                mol.particles.add(s.particles[n])

            bonds = [x for x in mol.bonds]
            mol.bonds = ItemContainer()
            for n in bonds:
                mol.bonds.add(s.bonds[n])

            angles = [x for x in mol.angles]
            mol.angles = ItemContainer()
            for n in angles:
                mol.angles.add(s.angles[n])

            dihedrals = [x for x in mol.dihedrals]
            mol.dihedrals = ItemContainer()
            for n in dihedrals:
                mol.dihedrals.add(s.dihedrals[n])

            impropers = [x for x in mol.impropers]
            mol.impropers = ItemContainer()
            for n in impropers:
                mol.impropers.add(s.impropers[n])

            s.molecules.add(mol)

    for p in s.particles:
        if s.particle_types[p.type]:
            p.type = s.particle_types[p.type]

        if s.molecules[p.molecule]:
            p.molecule = s.molecules[p.molecule]

        bonds = [x for x in p.bonds]
        p.bonds = ItemContainer()
        for n in bonds:
            p.bonds.add(s.bonds[n])

        angles = [x for x in p.angles]
        for n in angles:
            p.angles.add(s.angles[n])

        dihedrals = [x for x in p.dihedrals]
        for n in dihedrals:
            p.dihedrals.add(s.dihedrals[n])

        impropers = [x for x in p.impropers]
        for n in impropers:
            p.impropers.add(s.impropers[n])

    for b in s.bonds:
        if s.bond_types[b.type]:
            b.type = s.bond_types[b.type]
        b.a = s.particles[b.a]
        b.b = s.particles[b.b]

    for a in s.angles:
        if s.angle_types[a.type]:
            a.type = s.angle_types[a.type]
        a.a = s.particles[a.a]
        a.b = s.particles[a.b]
        a.c = s.particles[a.c]

    for d in s.dihedrals:
        if s.dihedral_types[d.type]:
            d.type = s.dihedral_types[d.type]
        d.a = s.particles[d.a]
        d.b = s.particles[d.b]
        d.c = s.particles[d.c]
        d.d = s.particles[d.d]

    for i in s.impropers:
        if s.improper_types[i.type]:
            i.type = s.improper_types[i.type]
        i.a = s.particles[i.a]
        i.b = s.particles[i.b]
        i.c = s.particles[i.c]
        i.d = s.particles[i.d]

    return s


def read_xyz(file_, **kwargs):
    """pysimm.system.read_xyz

    Interprets xyz file and creates pysimm.system.System object

    Args:
        file_: xyz file name
        quiet(optional): if False, print status

    Returns:
        pysimm.system.System object
    """
    quiet = kwargs.get('quiet')

    if os.path.isfile(file_):
        debug_print('reading file')
        f = file(file_)
    elif isinstance(file_, basestring):
        debug_print('reading string')
        f = StringIO(file_)

    s = System()
    nparticles = int(f.next().strip())
    name = f.next().strip()
    s.name = name
    for _ in range(nparticles):
        elem, x, y, z = f.next().split()
        x = float(x)
        y = float(y)
        z = float(z)
        s.particles.add(Particle(elem=elem, x=x, y=y, z=z))

    f.close()

    for p in s.particles:
        pt = s.particle_types.get(p.elem)
        if pt:
            p.type = pt[0]
        else:
            pt = ParticleType(elem=p.elem, name=p.elem)
            p.type = pt
            s.particle_types.add(pt)
    if not quiet:
        verbose_print('read %s particles' % s.particles.count)

    s.set_box(padding=0.5)

    return s


def read_lammps(data_file, **kwargs):
    """pysimm.system.read_lammps

    Interprets LAMMPS data file and creates pysimm.system.System object

    Args:
        data_file: LAMMPS data file name
        quiet(optional): if False, print status
        pair_style (optional): option to let user override
        bond_style (optional): option to let user override
        angle_style (optional): option to let user override
        dihedral_style (optional): option to let user override
        improper_style (optional): option to let user override
        set_types (optional): if True, objectify default=True
        name (optional): provide name for system

    Returns:
        pysimm.system.System object
    """
    pair_style = kwargs.get('pair_style')
    bond_style = kwargs.get('bond_style')
    angle_style = kwargs.get('angle_style')
    dihedral_style = kwargs.get('dihedral_style')
    improper_style = kwargs.get('improper_style')
    set_types = kwargs.get('set_types') or True
    name = kwargs.get('name')

    quiet = kwargs.get('quiet')

    if os.path.isfile(data_file):
        if not quiet:
            verbose_print('reading lammps data file "%s"' % data_file)
        f = file(data_file)
    elif isinstance(data_file, basestring):
        if not quiet:
            verbose_print('reading lammps data file from string')
        f = StringIO(data_file)
    else:
        raise PysimmError('pysimm.system.read_lammps requires either '
                          'file or string as first argument')

    if name:
        if not quiet:
            verbose_print('creating pysimm.system.System object with name %s'
                      % name)
        s = System(name=name)
    else:
        s = System(name=f.next().strip())

    nparticles = nparticle_types = nbonds = nbond_types = 0
    nangles = nangle_types = ndihedrals = ndihedral_types = 0
    nimpropers = nimproper_types = 0

    for line in f:
        line = line.split()
        if len(line) > 1 and line[1] == 'atoms':
            nparticles = int(line[0])
        elif len(line) > 1 and line[1] == 'atom':
            nparticle_types = int(line[0])
        elif len(line) > 1 and line[1] == 'bonds':
            nbonds = int(line[0])
        elif len(line) > 1 and line[1] == 'bond':
            nbond_types = int(line[0])
        elif len(line) > 1 and line[1] == 'angles':
            nangles = int(line[0])
        elif len(line) > 1 and line[1] == 'angle':
            nangle_types = int(line[0])
        elif len(line) > 1 and line[1] == 'dihedrals':
            ndihedrals = int(line[0])
        elif len(line) > 1 and line[1] == 'dihedral':
            ndihedral_types = int(line[0])
        elif len(line) > 1 and line[1] == 'impropers':
            nimpropers = int(line[0])
        elif len(line) > 1 and line[1] == 'improper':
            nimproper_types = int(line[0])
        elif len(line) > 3 and line[2] == 'xlo':
            s.dim.xlo = float(line[0])
            s.dim.xhi = float(line[1])
            s.dim.dx = s.dim.xhi - s.dim.xlo
        elif len(line) > 3 and line[2] == 'ylo':
            s.dim.ylo = float(line[0])
            s.dim.yhi = float(line[1])
            s.dim.dy = s.dim.yhi - s.dim.ylo
        elif len(line) > 3 and line[2] == 'zlo':
            s.dim.zlo = float(line[0])
            s.dim.zhi = float(line[1])
            s.dim.dz = s.dim.zhi - s.dim.zlo
        elif len(line) > 0 and line[0] == 'Masses':
            f.next()
            for i in range(nparticle_types):
                line = f.next().split('#')
                if len(line) == 2:
                    line, name = line
                    name = ','.join(re.split(',|\s+', name.strip()))
                else:
                    line = line[0]
                    name = None
                line = line.strip().split()
                tag = int(line[0])
                if s.particle_types[tag]:
                    if name is not None:
                        s.particle_types[tag].name = name
                    s.particle_types[tag].mass = float(line[1])
                else:
                    s.particle_types.add(ParticleType(tag=tag, name=name,
                                                      mass=float(line[1])))
            if not quiet:
                verbose_print('read masses for %s ParticleTypes'
                              % s.particle_types.count)
        elif len(line) > 0 and line[0] == 'Pair':
            if '#' in line and not pair_style:
                line = ' '.join(line).split('#')
                pair_style = line[1].strip()
            f.next()
            for i in range(nparticle_types):
                line = f.next().split('#')
                if len(line) == 2:
                    line, name = line
                    name = ','.join(re.split(',|\s+', name.strip()))
                else:
                    line = line[0]
                    name = None
                line = line.strip().split()
                tag = int(line[0])
                if pair_style and (pair_style.lower().startswith('lj') or
                                   pair_style.lower().startswith('class2')):
                    if s.particle_types[tag]:
                        pt = s.particle_types[tag]
                        if name is not None and pt.name is None:
                            pt.name = name
                        pt.epsilon = float(line[1])
                        pt.sigma = float(line[2])
                    else:
                        pt = ParticleType(tag=tag, name=name,
                                          epsilon=float(line[1]),
                                          sigma=float(line[2]))
                        s.particle_types.add(pt)
                elif pair_style and pair_style.lower().startswith('buck'):
                    if s.particle_types[tag]:
                        pt = s.particle_types[tag]
                        if name is not None and pt.name is None:
                            pt.name = name
                        pt.a = float(line[1])
                        pt.rho = float(line[2])
                        pt.c = float(line[3])
                    else:
                        pt = ParticleType(tag=tag, name=name, a=float(line[1]),
                                          rho=float(line[2]), c=float(line[3]))
                        s.particle_types.add(pt)
                elif not pair_style:
                    if not quiet and i == 0:
                        warning_print('pair_style not explicitly provided - '
                                      'guessing based on number of parameters '
                                      '(2=lj 3=buck)')
                    if len(line) == 3:
                        pair_style = 'lj'
                        if s.particle_types[tag]:
                            pt = s.particle_types[tag]
                            if name is not None and pt.name is None:
                                pt.name = name
                            pt.epsilon = float(line[1])
                            pt.sigma = float(line[2])
                        else:
                            pt = ParticleType(tag=tag, name=name,
                                              epsilon=float(line[1]),
                                              sigma=float(line[2]))
                            s.particle_types.add(pt)
                    elif len(line) == 4:
                        pair_style = 'buckingham'
                        if s.particle_types[tag]:
                            pt = s.particle_types[tag]
                            if name is not None and pt.name is None:
                                pt.name = name
                            pt.a = float(line[1])
                            pt.rho = float(line[2])
                            pt.c = float(line[3])
                        else:
                            pt = ParticleType(tag=tag, name=name, a=float(line[1]),
                                              rho=float(line[2]), c=float(line[3]))
                            s.particle_types.add(pt)
            if not quiet and pair_style:
                verbose_print('read "%s" nonbonded parameters '
                              'for %s ParticleTypes'
                              % (pair_style, nparticle_types))
            elif not quiet and not pair_style:
                verbose_print('cannot determine pair_style - '
                              'skipping nonbonded parameters')
        elif len(line) > 0 and line[0] == 'Bond':
            if '#' in line and not bond_style:
                line = ' '.join(line).split('#')
                bond_style = line[1].strip()
            f.next()
            for i in range(nbond_types):
                line = f.next().split('#')
                if len(line) == 2:
                    line, name = line
                    name = ','.join(re.split(',|\s+', name.strip()))
                else:
                    line = line[0]
                    name = None
                line = line.strip().split()
                tag = int(line[0])
                if bond_style and bond_style.lower().startswith('harm'):
                    s.bond_types.add(BondType(tag=tag, name=name,
                                              k=float(line[1]),
                                              r0=float(line[2])))
                elif bond_style and bond_style.lower().startswith('class2'):
                    s.bond_types.add(BondType(tag=tag, name=name,
                                              r0=float(line[1]),
                                              k2=float(line[2]),
                                              k3=float(line[3]),
                                              k4=float(line[4])))
                elif not bond_style:
                    if not quiet and i == 0:
                        warning_print('bond_style currently unknown - '
                                      'guessing based on number of parameters '
                                      '(2=harmonic 4=class2)')
                    if len(line) == 3:
                        bond_style = 'harmonic'
                        s.bond_types.add(BondType(tag=tag, name=name,
                                                  k=float(line[1]),
                                                  r0=float(line[2])))
                    elif len(line) == 5:
                        bond_style = 'class2'
                        s.bond_types.add(BondType(tag=tag, name=name,
                                                  r0=float(line[1]),
                                                  k2=float(line[2]),
                                                  k3=float(line[3]),
                                                  k4=float(line[4])))
            if not quiet and bond_style:
                verbose_print('read "%s" bond parameters '
                              'for %s BondTypes'
                              % (bond_style, nbond_types))
            elif not quiet and not bond_style:
                verbose_print('cannot determine bond_style - '
                              'skipping bond parameters')
        elif len(line) > 0 and line[0] == 'Angle':
            if '#' in line and not angle_style:
                line = ' '.join(line).split('#')
                angle_style = line[1].strip()
            f.next()
            for i in range(nangle_types):
                line = f.next().split('#')
                if len(line) == 2:
                    line, name = line
                    name = ','.join(re.split(',|\s+', name.strip()))
                else:
                    line = line[0]
                    name = None
                line = line.strip().split()
                tag = int(line[0])
                if angle_style and angle_style.lower().startswith('harm'):
                    s.angle_types.add(AngleType(tag=tag, name=name,
                                                k=float(line[1]),
                                                theta0=float(line[2])))
                elif angle_style and angle_style.lower().startswith('class2'):
                    s.angle_types.add(AngleType(tag=tag, name=name,
                                                theta0=float(line[1]),
                                                k2=float(line[2]),
                                                k3=float(line[3]),
                                                k4=float(line[4])))
                elif not angle_style:
                    if not quiet and i == 0:
                        warning_print('angle_style currently unknown - '
                                      'guessing based on number of parameters '
                                      '(2=harmonic 4=class2)')
                    if len(line) == 3:
                        angle_style = 'harmonic'
                        s.angle_types.add(AngleType(tag=tag, name=name,
                                                    k=float(line[1]),
                                                    theta0=float(line[2])))
                    elif len(line) == 5:
                        angle_style = 'class2'
                        s.angle_types.add(AngleType(tag=tag, name=name,
                                                    theta0=float(line[1]),
                                                    k2=float(line[2]),
                                                    k3=float(line[3]),
                                                    k4=float(line[4])))
            if not quiet and angle_style:
                    verbose_print('read "%s" angle parameters '
                                  'for %s AngleTypes'
                                  % (angle_style, nangle_types))
            elif not quiet and not angle_style:
                verbose_print('cannot determine angle_style - '
                              'skipping angle parameters')
        elif len(line) > 0 and line[0] == 'BondBond':
            f.next()
            for i in range(nangle_types):
                line = f.next().strip().split()
                tag = int(line[0])
                s.angle_types[tag].m = float(line[1])
                s.angle_types[tag].r1 = float(line[2])
                s.angle_types[tag].r2 = float(line[3])
            if not quiet and angle_style:
                verbose_print('read "%s" angle (bond-bond) '
                              'parameters for %s AngleTypes'
                              % (angle_style, nangle_types))
        elif len(line) > 0 and line[0] == 'BondAngle':
            f.next()
            for i in range(nangle_types):
                line = f.next().strip().split()
                tag = int(line[0])
                s.angle_types[tag].n1 = float(line[1])
                s.angle_types[tag].n2 = float(line[2])
            if not quiet and angle_style:
                verbose_print('read "%s" angle (bond-angle) '
                              'parameters for %s AngleTypes'
                              % (angle_style, nangle_types))
        elif len(line) > 0 and line[0] == 'Dihedral':
            if '#' in line and not dihedral_style:
                line = ' '.join(line).split('#')
                dihedral_style = line[1].strip()
            f.next()
            for i in range(ndihedral_types):
                line = f.next().split('#')
                if len(line) == 2:
                    line, name = line
                    name = ','.join(re.split(',|\s+', name.strip()))
                else:
                    line = line[0]
                    name = None
                line = line.strip().split()
                tag = int(line[0])
                if dihedral_style and dihedral_style.lower().startswith('harm'):
                    s.dihedral_types.add(DihedralType(tag=tag, name=name,
                                                      k=float(line[1]),
                                                      d=int(line[2]),
                                                      n=int(line[3])))
                elif (dihedral_style and
                        dihedral_style.lower().startswith('class2')):
                    s.dihedral_types.add(DihedralType(tag=tag, name=name,
                                                      k1=float(line[1]),
                                                      phi1=float(line[2]),
                                                      k2=float(line[3]),
                                                      phi2=float(line[4]),
                                                      k3=float(line[5]),
                                                      phi3=float(line[6])))
                elif not dihedral_style:
                    if not quiet and i == 0:
                        warning_print('dihedral_style currently unknown - '
                                      'guessing based on number of parameters '
                                      '(3=harmonic 6=class2)')
                    if len(line) == 4:
                        dihedral_style = 'harmonic'
                        s.dihedral_types.add(DihedralType(tag=tag, name=name,
                                                          k=float(line[1]),
                                                          d=int(line[2]),
                                                          n=int(line[3])))
                    elif len(line) == 7:
                        dihedral_style = 'class2'
                        s.dihedral_types.add(DihedralType(tag=tag, name=name,
                                                          k1=float(line[1]),
                                                          phi1=float(line[2]),
                                                          k2=float(line[3]),
                                                          phi2=float(line[4]),
                                                          k3=float(line[5]),
                                                          phi3=float(line[6])))
            if not quiet and dihedral_style:
                verbose_print('read "%s" dihedral parameters '
                              'for %s DihedralTypes'
                              % (dihedral_style, ndihedral_types))
            elif not quiet and not dihedral_style:
                verbose_print('cannot determine dihedral_style - '
                              'skipping bond parameters')
        elif len(line) > 0 and line[0] == 'MiddleBondTorsion':
            f.next()
            for i in range(ndihedral_types):
                line = f.next().strip().split()
                tag = int(line[0])
                s.dihedral_types[tag].a1 = float(line[1])
                s.dihedral_types[tag].a2 = float(line[2])
                s.dihedral_types[tag].a3 = float(line[3])
                s.dihedral_types[tag].r2 = float(line[4])
            if not quiet and dihedral_style:
                verbose_print('read "%s" dihedral '
                              '(middle-bond-torsion parameters for '
                              '%s DihedralTypes'
                              % (dihedral_style, ndihedral_types))
        elif len(line) > 0 and line[0] == 'EndBondTorsion':
            f.next()
            for i in range(ndihedral_types):
                line = f.next().strip().split()
                tag = int(line[0])
                s.dihedral_types[tag].b1 = float(line[1])
                s.dihedral_types[tag].b2 = float(line[2])
                s.dihedral_types[tag].b3 = float(line[3])
                s.dihedral_types[tag].c1 = float(line[4])
                s.dihedral_types[tag].c2 = float(line[5])
                s.dihedral_types[tag].c3 = float(line[6])
                s.dihedral_types[tag].r1 = float(line[7])
                s.dihedral_types[tag].r3 = float(line[8])
            if not quiet and dihedral_style:
                verbose_print('read "%s" dihedral '
                              '(end-bond-torsion parameters for '
                              '%s DihedralTypes'
                              % (dihedral_style, ndihedral_types))
        elif len(line) > 0 and line[0] == 'AngleTorsion':
            f.next()
            for i in range(ndihedral_types):
                line = f.next().strip().split()
                tag = int(line[0])
                s.dihedral_types[tag].d1 = float(line[1])
                s.dihedral_types[tag].d2 = float(line[2])
                s.dihedral_types[tag].d3 = float(line[3])
                s.dihedral_types[tag].e1 = float(line[4])
                s.dihedral_types[tag].e2 = float(line[5])
                s.dihedral_types[tag].e3 = float(line[6])
                s.dihedral_types[tag].theta1 = float(line[7])
                s.dihedral_types[tag].theta2 = float(line[8])
            if not quiet and dihedral_style:
                    verbose_print('read "%s" dihedral '
                                  '(angle-torsion parameters for '
                                  '%s DihedralTypes'
                                  % (dihedral_style, ndihedral_types))
        elif len(line) > 0 and line[0] == 'AngleAngleTorsion':
            f.next()
            for i in range(ndihedral_types):
                line = f.next().strip().split()
                tag = int(line[0])
                s.dihedral_types[tag].m = float(line[1])
            if not quiet and dihedral_style:
                verbose_print('read "%s" dihedral '
                              '(angle-angle-torsion parameters for '
                              '%s DihedralTypes'
                              % (dihedral_style, ndihedral_types))
        elif len(line) > 0 and line[0] == 'BondBond13':
            f.next()
            for i in range(ndihedral_types):
                line = f.next().strip().split()
                tag = int(line[0])
                s.dihedral_types[tag].n_class2 = float(line[1])
            if not quiet and dihedral_style:
                verbose_print('read "%s" dihedral '
                              '(bond-bond-1-3 parameters for '
                              '%s DihedralTypes'
                              % (dihedral_style, ndihedral_types))
        elif len(line) > 0 and line[0] == 'Improper':
            if '#' in line and not improper_style:
                line = ' '.join(line).split('#')
                improper_style = line[1].strip()
            f.next()
            for i in range(nimproper_types):
                line = f.next().split('#')
                if len(line) == 2:
                    line, name = line
                    name = ','.join(re.split(',|\s+', name.strip()))
                else:
                    line = line[0]
                    name = None
                line = line.strip().split()
                tag = int(line[0])
                if improper_style and improper_style.lower().startswith('harm'):
                    s.improper_types.add(ImproperType(tag=tag, name=name,
                                                      k=float(line[1]),
                                                      x0=float(line[2])))
                elif (improper_style and
                      improper_style.lower().startswith('class2')):
                    s.improper_types.add(ImproperType(tag=tag, name=name,
                                                      k=float(line[1]),
                                                      x0=float(line[2])))
                elif not improper_style:
                    if not quiet and i == 0:
                            warning_print('cannot guess improper_style '
                                          'from number of parameters - '
                                          'will try to determine style later '
                                          'based on other types')
                    s.improper_types.add(ImproperType(tag=tag, name=name,
                                                      k=float(line[1]),
                                                      x0=float(line[2])))

            if not quiet and improper_style:
                verbose_print('read "%s" improper parameters '
                              'for %s ImproperTypes'
                              % (improper_style, nimproper_types))
        elif len(line) > 0 and line[0] == 'AngleAngle':
            f.next()
            for i in range(nimproper_types):
                line = f.next().strip().split()
                tag = int(line[0])
                s.improper_types[tag].m1 = float(line[1])
                s.improper_types[tag].m2 = float(line[2])
                s.improper_types[tag].m3 = float(line[3])
                s.improper_types[tag].theta1 = float(line[4])
                s.improper_types[tag].theta2 = float(line[5])
                s.improper_types[tag].theta3 = float(line[6])
            if not quiet and improper_style:
                verbose_print('read "%s" improper '
                              '(angle-angle parameters for '
                              '%s ImproperTypes'
                              % (improper_style, nimproper_types))
        elif len(line) > 0 and line[0] == 'Atoms':
            f.next()
            for i in range(nparticles):
                line = f.next().strip().split()
                tag = int(line[0])
                if s.particles[tag]:
                    p = s.particles[tag]
                    d_ = {'molecule': int(line[1]), 'type': int(line[2]),
                          'charge': float(line[3]), 'x': float(line[4]),
                          'y': float(line[5]), 'z': float(line[6])}
                    p.set(**d_)
                else:
                    p = Particle(tag=tag, molecule=int(line[1]),
                                 type=int(line[2]),
                                 charge=float(line[3]), x=float(line[4]),
                                 y=float(line[5]), z=float(line[6]),
                                 vx=0., vy=0., vz=0.)
                    s.particles.add(p)
                if s.dim.check():
                    p.frac_x = p.x / s.dim.dx
                    p.frac_y = p.y / s.dim.dy
                    p.frac_z = p.z / s.dim.dz
            if not quiet:
                verbose_print('read %s particles' % nparticles)
        elif len(line) > 0 and line[0] == 'Velocities':
            f.next()
            for i in range(nparticles):
                line = f.next().strip().split()
                tag = int(line[0])
                if s.particles[tag]:
                    p = s.particles[tag]
                    d_ = {'vx': float(line[1]), 'vy': float(line[2]),
                          'vz': float(line[3])}
                    p.set(**d_)
                else:
                    p = Particle(tag=tag, vx=float(line[1]), vy=float(line[2]),
                                 vz=float(line[3]))
                    s.particles.add(p)
            if not quiet:
                verbose_print('read velocities for %s particles' % nparticles)
        elif len(line) > 0 and line[0] == 'Bonds':
            f.next()
            for i in range(nbonds):
                line = f.next().strip().split()
                tag = int(line[0])
                b = Bond(tag=tag, type=int(line[1]),
                         a=int(line[2]), b=int(line[3]))
                s.bonds.add(b)
            if not quiet:
                verbose_print('read %s bonds' % nbonds)
        elif len(line) > 0 and line[0] == 'Angles':
            f.next()
            for i in range(nangles):
                line = f.next().strip().split()
                tag = int(line[0])
                a = Angle(tag=tag, type=int(line[1]),
                          a=int(line[2]), b=int(line[3]), c=int(line[4]))
                s.angles.add(a)
            if not quiet:
                verbose_print('read %s angles' % nangles)
        elif len(line) > 0 and line[0] == 'Dihedrals':
            f.next()
            for i in range(ndihedrals):
                line = f.next().strip().split()
                tag = int(line[0])
                d = Dihedral(tag=tag, type=int(line[1]),
                             a=int(line[2]), b=int(line[3]),
                             c=int(line[4]), d=int(line[5]))
                s.dihedrals.add(d)
            if not quiet:
                verbose_print('read %s dihedrals' % ndihedrals)
        elif len(line) > 0 and line[0] == 'Impropers':
            f.next()
            for i in range(nimpropers):
                line = f.next().strip().split()
                tag = int(line[0])
                if (s.ff_class == '2' or improper_style == 'class2' or (s.improper_types[1] and s.improper_types[1].m1
                        is not None)):
                    s.impropers.add(Improper(tag=tag, type=int(line[1]),
                                             a=int(line[3]), b=int(line[2]),
                                             c=int(line[4]), d=int(line[5])))
                else:
                    s.impropers.add(Improper(tag=tag, type=int(line[1]),
                                             a=int(line[2]), b=int(line[3]),
                                             c=int(line[4]), d=int(line[5])))
            if not quiet:
                verbose_print('read %s impropers' % nimpropers)
    f.close()

    s.pair_style = pair_style
    s.bond_style = bond_style
    s.angle_style = angle_style
    s.dihedral_style = dihedral_style
    if improper_style:
        s.improper_style = improper_style
    elif not improper_style and s.impropers.count > 1:
        if not quiet:
            verbose_print('improper style not set explicitly '
                          'but impropers exist in system, guessing style '
                          'based on other forcefield styles...')
        if (s.bond_style.startswith('harm') or
                s.angle_style.startswith('harm') or
                s.dihedral_style.startswith('harm')):
            improper_style = 'harmonic'
            s.improper_style = 'harmonic'
        elif (s.bond_style.startswith('class2') or
                s.angle_style.startswith('class2') or
                s.dihedral_style.startswith('class2')):
            improper_style = 'class2'
            s.improper_style = 'class2'
        if s.improper_style:
            if not quiet:
                verbose_print('setting improper style to "%s", '
                              'if this is incorrect try explicitly setting '
                              'improper_style as argument in '
                              'pysimm.system.read_lammps' % improper_style)
        else:
            if not quiet:
                error_print('still cannot determine improper style...')

    if pair_style and pair_style.startswith('lj'):
        if ((s.bond_style and s.bond_style.startswith('class2')) or
                (s.angle_style and s.angle_style.startswith('class2')) or
                (s.dihedral_style and s.dihedral_style.startswith('class2'))):
            s.pair_style = 'class2'

    styles = [s.pair_style, s.bond_style, s.angle_style, s.dihedral_style,
              s.improper_style]
    if 'class2' in styles:
        s.ff_class = '2'
    else:
        s.ff_class = '1'

    if 'harmonic' in styles and 'class2' in styles:
        if not quiet:
            warning_print('it appears there is a mixture of class1 and class2 '
                          'forcefield styles in your system...this is usually '
                          'unadvised')

    if set_types:
        s.objectify()

    for pt in s.particle_types:
        if pt.name and pt.name.find('@') >= 0:
            if pt.name.split('@')[-1][0].upper() in ['H', 'C', 'N', 'O', 'F', 'S']:
                pt.elem = pt.name.split('@')[-1][0].upper()
        if pt.name and pt.name[0] == 'L' and pt.name[1] != 'i':
            pt.elem = pt.name[1].upper()
        elif pt.name:
            if pt.name[1:3] == 'Na':
                pt.elem = 'Na'
            if pt.name[0].upper() in ['H', 'C', 'N', 'O', 'F', 'S']:
                pt.elem = pt.name[0].upper()

    for p in s.particles:
        if isinstance(p.type, ParticleType) and p.type.name and p.type.name.find('@') >= 0:
            if p.type.name[0].upper() == 'H':
                p.linker = 'head'
            elif p.type.name[0].upper() == 'T':
                p.linker = 'tail'
            elif p.type.name[0].upper() == 'L':
                p.linker = True

    s.set_cog()
    s.set_mass()
    s.set_volume()
    s.set_density()
    s.set_velocity()

    return s


def read_hoomd(data_file, **kwargs):
    """pysimm.system.read_hoomd

    Interprets hoomd data file and creates pysimm.system.System object

    Args:
        data_file: hoomd data file name
        f (optional): Forcefield object to get data from
        d_unit (optional): allow user to override distance unit used in data file default='nm'

    Returns:
        pysimm.system.System object
    """
    f = kwargs.get('forcefield')
    string = kwargs.get('string')
    d_unit = kwargs.get('d_unit') or 'nm'

    s = System()
    if string:
        root = Et.fromstring(string)
    else:
        tree = Et.parse(data_file)
        root = tree.getroot()
    config = root.find('configuration')
    box = config.find('box')
    dims = box.attrib
    s.dim.xlo = -1*float(dims.get('lx'))/2
    s.dim.xhi = float(dims.get('lx'))/2
    s.dim.ylo = -1*float(dims.get('ly'))/2
    s.dim.yhi = float(dims.get('ly'))/2
    s.dim.zlo = -1*float(dims.get('lz'))/2
    s.dim.zhi = float(dims.get('lz'))/2
    if d_unit == 'nm':
        s.dim.xlo *= 10
        s.dim.xhi *= 10
        s.dim.ylo *= 10
        s.dim.yhi *= 10
        s.dim.zlo *= 10
        s.dim.zhi *= 10
    position = config.find('position').text.split()
    mass = config.find('mass').text.split() if config.find('mass') is not None else None
    charge = config.find('charge').text.split() if config.find('charge') is not None else None
    types = config.find('type').text.split()
    bond = config.find('bond').text.split() if config.find('bond') is not None else None
    angle = config.find('angle').text.split() if config.find('angle') is not None else None
    dihedral = config.find('dihedral').text.split() if config.find('dihedral') is not None else None
    improper = config.find('improper').text.split() if config.find('improper') is not None else None

    nparticles = len(types)
    nbonds = len(bond)/3 if bond else 0
    nangles = len(angle)/4 if angle else 0
    ndihedrals = len(dihedral)/5 if dihedral else 0
    nimpropers = len(improper)/5 if improper else 0

    ptypes = set()
    for i in range(nparticles):
        ptypes.add(types[i])
        x = float(position[3*i])
        y = float(position[3*i+1])
        z = float(position[3*i+2])
        if d_unit == 'nm':
            x *= 10
            y *= 10
            z *= 10
        p = Particle(tag=i+1, type=types[i], molecule=1, x=x, y=y, z=z, mass=mass[i])
        s.particles.add(p)
        if charge:
            p.charge = charge[i]
        else:
            p.charge = 0.
    ptypes = list(ptypes)

    btypes = set()
    for i in range(nbonds):
        btypes.add(bond[3*i])
        s.bonds.add(Bond(tag=i+1, type=bond[3*i], a=int(bond[3*i+1])+1, b=int(bond[3*i+2])+1))
    btypes = list(btypes)

    atypes = set()
    for i in range(nangles):
        atypes.add(angle[4*i])
        s.angles.add(Angle(tag=i+1, type=angle[4*i], a=int(angle[4*i+1])+1, b=int(angle[4*i+2])+1,
                           c=int(angle[4*i+3])+1))
    atypes = list(atypes)

    dtypes = set()
    for i in range(ndihedrals):
        dtypes.add(dihedral[5*i])
        s.dihedrals.add(Dihedral(tag=i+1, type=dihedral[5*i], a=int(dihedral[5*i+1])+1, b=int(dihedral[5*i+2])+1,
                                 c=int(dihedral[5*i+3])+1, d=int(dihedral[5*i+4])+1))
    dtypes = list(dtypes)

    itypes = set()
    for i in range(nimpropers):
        itypes.add(improper[5*i])
        s.impropers.add(Improper(tag=i+1, type=improper[5*i], a=int(improper[5*i+1])+1, b=int(improper[5*i+2])+1,
                                 c=int(improper[5*i+3])+1, d=int(improper[5*i+4])+1))
    itypes = list(itypes)

    nparticle_types = len(ptypes)
    nbond_types = len(btypes)
    nangle_types = len(atypes)
    ndihedral_types = len(dtypes)
    nimproper_types = len(itypes)

    for i in range(nparticle_types):
        s.particle_types.add(ParticleType(tag=i+1, name=ptypes[i]))

    for i in range(nbond_types):
        s.bond_types.add(BondType(tag=i+1, name=btypes[i]))

    for i in range(nangle_types):
        s.angle_types.add(AngleType(tag=i+1, name=atypes[i]))

    for i in range(ndihedral_types):
        s.dihedral_types.add(DihedralType(tag=i+1, name=dtypes[i]))

    for i in range(nimproper_types):
        s.improper_types.add(ImproperType(tag=i+1, name=itypes[i]))

    if f:

        for pt in s.particle_types:
            if f.particle_types.get(pt.name):
                tag = pt.tag
                s.particle_types.remove(tag, update=False)
                pt = f.particle_types.get(pt.name)[0].copy()
                pt.tag = tag
                s.particle_types.add(pt)

        for bt in s.bond_types:
            if f.bond_types.get(bt.name):
                tag = bt.tag
                s.bond_types.remove(tag, update=False)
                bt = f.bond_types.get(bt.name)[0].copy()
                bt.tag = tag
                s.bond_types.add(bt)

        for at in s.angle_types:
            if f.angle_types.get(at.name):
                tag = at.tag
                s.angle_types.remove(tag, update=False)
                at = f.angle_types.get(at.name)[0].copy()
                at.tag = tag
                s.angle_types.add(at)

        for dt in s.dihedral_types:
            if f.dihedral_types.get(dt.name):
                tag = dt.tag
                s.dihedral_types.remove(tag, update=False)
                dt = f.dihedral_types.get(dt.name)[0].copy()
                dt.tag = tag
                s.dihedral_types.add(dt)

        for it in s.improper_types:
            if f.improper_types.get(it.name):
                tag = it.tag
                s.improper_types.remove(tag, update=False)
                it = f.improper_types.get(it.name)[0].copy()
                it.tag = tag
                s.improper_types.add(it)

    else:
        warning_print('no forcefield supplied; types will be incomplete')

    for p in s.particles:
        for pt in s.particle_types:
            if p.type == pt.name:
                p.type = pt
                p.type.mass = float(p.mass)
                break

    for b in s.bonds:
        for bt in s.bond_types:
            if b.type == bt.name:
                b.type = bt
                break

    for a in s.angles:
        for at in s.angle_types:
            if a.type == at.name:
                a.type = at
                break

    for d in s.dihedrals:
        for dt in s.dihedral_types:
            if d.type == dt.name:
                d.type = dt
                break

    for i in s.impropers:
        for it in s.improper_types:
            if i.type == it.name:
                i.type = it
                break

    s.objectify()

    s.set_cog()
    s.set_mass()
    s.set_volume()
    s.set_density()
    s.set_velocity()

    return s


def read_pubchem_smiles(smiles, type_with=None):
    """pysimm.system.read_pubchem_smiles

    Interface with pubchem restful API to create molecular system from SMILES format

    Args:
        smiles: smiles formatted string of molecule
        type_with: pysimm.forcefield.Forcefield object to type with default=None

    Returns:
        pysimm.system.System object
    """

    print('making request to pubchem RESTful API...')

    req = ('https://pubchem.ncbi.nlm.nih.gov/'
           'rest/pug/compound/smiles/%s/SDF/?record_type=3d' % smiles)

    try:
        resp = urlopen(req)
        return read_mol(resp.read(), type_with=type_with)
    except HTTPError, URLError:
        print('Could not retrieve pubchem entry for smiles %s' % smiles)


def read_cml(cml_file, **kwargs):
    """pysimm.system.read_cml

    Interprets cml file and creates pysimm.system.System object

    Args:
        cml_file: cml file name
        linkers (optional): if True, use spinMultiplicity to determine linker default=None

    Returns:
        pysimm.system.System object
    """
    linkers = kwargs.get('linkers')

    if os.path.isfile(cml_file):
        debug_print('reading file')
        iter_parse = Et.iterparse(cml_file)
    elif isinstance(cml_file, basestring):
        debug_print('reading string')
        iter_parse = Et.iterparse(StringIO(cml_file))
    else:
        raise PysimmError('pysimm.system.read_cml requires a file as argument')

    for _, el in iter_parse:
        if '}' in el.tag:
            el.tag = el.tag.split('}', 1)[1]
    root = iter_parse.root

    s = System(name='read using pysimm.system.read_cml')

    particles = root.find('atomArray')
    bonds = root.find('bondArray')

    for p_ in particles:
        tag = int(p_.attrib['id'].replace('a', ''))
        elem = p_.attrib['elementType']
        x = float(p_.attrib['x3'])
        y = float(p_.attrib['y3'])
        z = float(p_.attrib['z3'])
        if linkers:
            linker = True if p_.attrib.get('spinMultiplicity') else None
        else:
            linker = None

        p = Particle(tag=tag, elem=elem, x=x, y=y, z=z, charge=0, molecule=1, linker=linker)
        s.particles.add(p)

    for b_ in bonds:
        a, b = b_.attrib['atomRefs2'].split()
        a = int(a.replace('a', ''))
        b = int(b.replace('a', ''))
        order = b_.attrib['order']
        if order == 'A':
            order = 4
        else:
            order = int(order)

        b = Bond(a=a, b=b, order=order)
        s.bonds.add(b)

    s.objectify()

    return s


def read_mol(mol_file, type_with=None, version='V2000'):
    """pysimm.system.read_mol

    Interprets mol file and creates pysimm.system.System object

    Args:
        mol_file: mol file name
        f (optional): Forcefield object to get data from
        version: version of mol file to expect default='V2000'

    Returns:
        pysimm.system.System object
    """
    if os.path.isfile(mol_file):
        debug_print('reading file')
        f = file(mol_file)
    elif isinstance(mol_file, basestring):
        debug_print('reading string')
        f = StringIO(mol_file)
    else:
        raise PysimmError('pysimm.system.read_mol requires either '
                          'file or string as argument')

    s = System(name='read using pysimm.system.read_mol')
    for n in range(3):
        f.next()

    line = f.next()

    nparticles = int(line.split()[0])
    nbonds = int(line.split()[1])
    if len(line.split()) >= 3:
        version = line.split()[-1]

    if version == 'V2000':
        for n in range(nparticles):
            line = f.next()
            x, y, z, elem, something, charge = line.split()[:6]
            p = Particle(x=float(x), y=float(y), z=float(z), molecule=1,
                         elem=elem, charge=float(charge))
            s.particles.add(p)
            if p.elem[0] == 'L':
                p.linker = True
                p.elem = p.elem[1:]
            elif p.charge == 5:
                p.linker = True
                p.charge = 0

        for n in range(nbonds):
            line = f.next()
            a, b, order = map(int, line.split()[:3])
            new_bond = s.bonds.add(Bond(a=a, b=b, order=order))

    elif version == 'V3000':
        f.next()
        line = f.next()
        nparticles = int(line.split()[3])
        nbonds = int(line.split()[4])
        f.next()

        for n in range(nparticles):
            line = f.next()
            id_, elem, x, y, z, charge = line.split()[2:8]
            p = Particle(x=float(x), y=float(y), z=float(z), molecule=1,
                         elem=elem, charge=float(charge))
            s.particles.add(p)

        f.next()
        f.next()

        for n in range(nbonds):
            line = f.next()
            id_, order, a, b = map(int, line.split()[2:6])
            s.bonds.add(Bond(a=a, b=b, order=order))

    s.objectify()

    if type_with:
        try:
            s.apply_forcefield(type_with)
        except Exception:
            print('forcefield typing with forcefield %s unsuccessful'
                  % type_with.ff_name)

    return s


def read_pdb(pdb_file, guess_bonds=False):
    """pysimm.system.read_pdb

    Interprets pdb file and creates pysimm.system.System object

    Args:
        pdb_file: pdb file name
        guess_bonds: if True, guess bonds default=False ** probably buggy ***

    Returns:
        pysimm.system.System object
    """
    if os.path.isfile(pdb_file):
        debug_print('reading file')
        f = file(pdb_file)
    elif isinstance(pdb_file, basestring):
        debug_print('reading string')
        f = StringIO(pdb_file)
    else:
        raise PysimmError('pysimm.system.read_pdb requires either '
                          'file or string as argument')

    s = System(name='read using pysimm.system.read_pdb')

    for line in f:
        if line.startswith('ATOM'):
            tag = int(line[6:11].strip())
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chainid = line[21]
            resid = line[22:26].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            elem = line[76:78].strip()
            p = Particle(tag=tag, name=name, resname=resname, chainid=chainid, resid=resid, x=x, y=y, z=z, elem=elem)
            if not s.particles[tag]:
                s.particles.add(p)


    f.close()

    if guess_bonds:
        s.guess_bonds()

    return s


def compare(s1, s2):
    print('Particle Types:\n')
    for pt in s1.particle_types:
        s2_pt = s2.particle_types.get(pt.name)
        if s2_pt and len(s2_pt) == 1:
            s2_pt = s2_pt[0]
            print('%s\n%s\n' % (vars(pt), vars(s2_pt)))

    print('\n\nBond Types:\n')
    for bt in s1.bond_types:
        s2_bt = s2.bond_types.get(bt.name)
        if s2_bt and len(s2_bt) == 1:
            s2_bt = s2_bt[0]
            print('%s\n%s\n' % (vars(bt), vars(s2_bt)))

    print('\n\nAngle Types:\n')
    for at in s1.angle_types:
        s2_at = s2.angle_types.get(at.name)
        if s2_at and len(s2_at) == 1:
            s2_at = s2_at[0]
            print('%s\n%s\n' % (vars(at), vars(s2_at)))

    print('\n\nDihedral Types:\n')
    for dt in s1.dihedral_types:
        s2_dt = s2.dihedral_types.get(dt.name)
        if s2_dt and len(s2_dt) == 1:
            s2_dt = s2_dt[0]
            print('%s\n%s\n' % (vars(dt), vars(s2_dt)))

    print('\n\nImproper Types:\n')
    for it in s1.improper_types:
        s2_it = s2.improper_types.get(it.name)
        if s2_it and len(s2_it) == 1:
            s2_it = s2_it[0]
            print('%s\n%s\n' % (vars(it), vars(s2_it)))


def get_types(*arg, **kwargs):
    """pysimm.system.get_types

    Get unique type names from list of systems

    Args:
        write (optional): if True, write types dictionary to filename

    Returns:
        (ptypes, btypes, atypes, dtypes, itypes)
        *** for use with update_types ***
    """
    write = kwargs.get('write')

    ptypes = ItemContainer()
    btypes = ItemContainer()
    atypes = ItemContainer()
    dtypes = ItemContainer()
    itypes = ItemContainer()

    for s in arg:
        for t in s.particle_types:
            if t.name and t.name not in [x.name for x in ptypes]:
                ptypes.add(t.copy())
        for t in s.bond_types:
            if t.name and t.name not in [x.name for x in btypes]:
                btypes.add(t.copy())
        for t in s.angle_types:
            if t.name and t.name not in [x.name for x in atypes]:
                atypes.add(t.copy())
        for t in s.dihedral_types:
            if t.name and t.name not in [x.name for x in dtypes]:
                dtypes.add(t.copy())
        for t in s.improper_types:
            if t.name and t.name not in [x.name for x in itypes]:
                itypes.add(t.copy())

    if write:
        t_file = open('types.txt', 'w+')
        if ptypes.count > 0:
            t_file.write('atom types\n')
        for t in ptypes:
            t_file.write('%s %s\n' % (t.tag, t.name))
        if btypes.count > 0:
            t_file.write('\nbond types\n')
        for t in btypes:
            t_file.write('%s %s\n' % (t.tag, t.name))
        if atypes.count > 0:
            t_file.write('\nangle types\n')
        for t in atypes:
            t_file.write('%s %s\n' % (t.tag, t.name))
        if dtypes.count > 0:
            t_file.write('\ndihedral types\n')
        for t in dtypes:
            t_file.write('%s %s\n' % (t.tag, t.name))
        if itypes.count > 0:
            t_file.write('\nimproper types\n')
        for t in itypes:
            t_file.write('%s %s\n' % (t.tag, t.name))
        t_file.close()

    return ptypes, btypes, atypes, dtypes, itypes


def distance_to_origin(p):
    return sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2))


def cluster(s, cutoff=35):
    s.clusters = ItemContainer()
    for p in s.particles:
        if p.cluster:
            for p0 in s.particles:
                if p0 is not p and calc.pbc_distance(p, p0) < cutoff:
                    p.cluster.add(p0)
        else:
            p.cluster = ItemContainer()
            s.clusters.add(p.cluster)
            for p0 in s.particles:
                if p0 is not p and calc.pbc_distance(p, p0) < cutoff:
                    p.cluster.add(p0)


def replicate(ref, nrep, s_=None, density=0.3, rand=True, print_insertions=True):
    """pysimm.system.replicate

    Replicates list of System objects into new (or exisintg) System.
    Can be random insertion.

    Args:
        ref: reference System(s) (this can be a list)
        nrep: number of insertions to perform (can be list but must match length of ref)
        s_: System into which insertions will be performed default=None
        density: density of new System default=0.3 (set to None to not change box)
        rand: if True, random insertion is performed
        print_insertions: if True, update screen with number of insertions
    """
    if not isinstance(ref, list):
        ref = [ref]
    if not isinstance(nrep, list):
        nrep = [nrep]
    assert len(ref) == len(nrep)

    if s_ is None:
        s_ = System()
    s_.ff_class = ref[0].ff_class
    s_.pair_style = ref[0].pair_style

    for r in ref:
        r.set_mass()
        r.center_at_origin()
        r.r = 0
        for p in r.particles:
            r.r = max(r.r, distance_to_origin(p))
        s_.molecule_types.add(r)

    mass = 0
    for i, r in enumerate(ref):
        mass += r.mass * nrep[i]
    mass /= 6.02e23

    if density:
        volume = float(mass) / density
        boxl = pow(volume, 1 / 3.) * 1e8
        s_.dim.xlo = -1. * boxl / 2.
        s_.dim.xhi = boxl / 2.
        s_.dim.dx = s_.dim.xhi - s_.dim.xlo
        s_.dim.ylo = -1. * boxl / 2.
        s_.dim.yhi = boxl / 2.
        s_.dim.dy = s_.dim.yhi - s_.dim.ylo
        s_.dim.zlo = -1. * boxl / 2.
        s_.dim.zhi = boxl / 2.
        s_.dim.dz = s_.dim.zhi - s_.dim.zlo

    num = 0
    for j, r in enumerate(ref):
        for n in range(nrep[j]):
            if rand:
                rotate_x = random() * 2 * pi
                rotate_y = random() * 2 * pi
                rotate_z = random() * 2 * pi
                dx = s_.dim.xhi - s_.dim.xlo
                dx = (-dx / 2. + r.r) + random() * (dx - 2 * r.r)
                dy = s_.dim.yhi - s_.dim.ylo
                dy = (-dy / 2. + r.r) + random() * (dy - 2 * r.r)
                dz = s_.dim.zhi - s_.dim.zlo
                dz = (-dz / 2. + r.r) + random() * (dz - 2 * r.r)
                r_ = r.copy(rotate_x=rotate_x, rotate_y=rotate_y,
                            rotate_z=rotate_z, dx=dx, dy=dy, dz=dz)
            else:
                r_ = r.copy()

            s_.add(r_, change_dim=False, update_properties=False)
            num += 1
            if print_insertions:
                verbose_print('Molecule %s inserted' % num)

    s_.set_density()
    s_.set_cog()
    s_.set_velocity()

    return s_
