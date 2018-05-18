# ******************************************************************************
# pysimm.forcefield.tip3p module
# ******************************************************************************
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

import os
from itertools import permutations, combinations

import gasteiger
from pysimm.system import Angle, Dihedral, Improper
from forcefield import Forcefield


class Tip3p(Forcefield):
    """pysimm.forcefield.Tip3p

    Forcefield object with typing rules for Tip3p model.
    By default reads data file in forcefields subdirectory.

    Attributes:
        ff_name: tip3p
        pair_style: lj
        ff_class: 1
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   os.pardir, os.pardir, 'dat', 'forcefields', 'tip3p.xml')
        Forcefield.__init__(self, db_file)
        self.name = 'tip3p'
        self.pair_style = 'lj'
        self.ff_class = '1'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Tip3p.assign_ptypes

        Tip3p specific particle typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds defined.
        *** use System.add_particle_bonding() to ensure this ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        all_types = set()
        s.pair_style = self.pair_style
        for p in s.particles:
            p.bonded_to = [x.a if p is x.b else x.b for x in p.bonds]
            p.bond_orders = [x.order for x in p.bonds]
            if None in p.bond_orders:
                error_print('error: bond orders are not set')
            p.bond_elements = [x.a.elem if p is x.b else x.b.elem for x in
                               p.bonds]
            p.nbonds = len(p.bond_elements)
        for p in s.particles:
            if p.elem == 'H':
                p.type_name = 'hw'
            elif p.elem == 'O':
                p.type_name = 'ow'
            else:
                print 'cant type particle %s' % p.tag
                return p

            type_ = self.particle_types.get(p.type_name)
            if not type_:
                print(p.tag, p.elem, p.type_name)

            all_types.add(self.particle_types.get(p.type_name)[0])

        for pt in all_types:
            pt = pt.copy()
            s.particle_types.add(pt)

        for p in s.particles:
            pt = s.particle_types.get(p.type_name)
            if pt:
                p.type = pt[0]

    def assign_btypes(self, s):
        """pysimm.forcefield.Tip3p.assign_btypes

        Tip3p specific bond typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have type and type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        all_types = set()
        for b in s.bonds:
            bt = self.bond_types.get('%s,%s' % (b.a.type.name, b.b.type.name))
            if bt:
                b.type_name = bt[0].name
            else:
                print ('couldnt type this bond %s,%s'
                       % (b.a.type.name, b.b.type.name))
                return b
            all_types.add(self.bond_types.get(b.type_name)[0])

        for bt in all_types:
            bt = bt.copy()
            s.bond_types.add(bt)

        for b in s.bonds:
            bt = s.bond_types.get(b.type_name)
            if bt:
                b.type = bt[0]

    def assign_atypes(self, s):
        """pysimm.forcefield.Tip3p.assign_atypes

        Tip3p specific angle typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds, type and type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        all_types = set()
        for p in s.particles:
            p.bonded_to = [x.a if p is x.b else x.b for x in p.bonds]
            for p1 in p.bonded_to:
                for p2 in p.bonded_to:
                    if p1 is not p2:
                        unique = True
                        for a in s.angles:
                            if ((a.a is p1 and a.b is p and a.c is p2) or
                                    (a.a is p2 and a.b is p and a.c is p1)):
                                unique = False
                        if unique:
                            at = self.angle_types.get('%s,%s,%s'
                                                      % (p1.type.name,
                                                         p.type.name,
                                                         p2.type.name))
                            if at:
                                s.angles.add(Angle(type_name=at[0].name,
                                                   a=p1, b=p, c=p2))
                                all_types.add(at[0])
                            else:
                                print ('I cant type this angle %s,%s,%s'
                                       % (p1.type.name,
                                          p.type.name,
                                          p2.type.name))

        for at in all_types:
            at = at.copy()
            s.angle_types.add(at)

        for a in s.angles:
            at = s.angle_types.get(a.type_name)
            if at:
                a.type = at[0]

    def assign_dtypes(self, s):
        """pysimm.forcefield.Tip3p.assign_dtypes

        Tip3p specific dihedral typing rules.
        There are none.

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        pass

    def assign_itypes(self, s):
        """pysimm.forcefield.Tip3p.assign_itypes

        Tip3p specific improper typing rules.
        There are none.

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        pass

    def assign_charges(self, s, charges='default'):
        """pysimm.forcefield.Tip3p.assign_charges

        Tip3p specific charge assignment.
        There are none.

        Args:
            s: :class:`~pysimm.system.System`
            charges: default

        Returns:
            None
        """
        if charges == 'gasteiger':
            print('adding gasteiger charges')
            gasteiger.set_charges(s)
        elif charges == 'default':
            print('adding default TIP3P charges')
            for p in s.particles:
                p.charge = 0
            for b in s.bonds:
                n1 = b.a.type.eq_bond or b.a.type.name
                n2 = b.b.type.eq_bond or b.b.type.name
                btype = self.bond_types.get('%s,%s' % (n1, n2))
                if btype:
                    btype = btype[0]
                    if btype.name == '%s,%s' % (n1, n2):
                        b.a.charge += float(btype.q1)
                        b.b.charge += float(btype.q2)
                    elif btype.name == '%s,%s' % (n2, n1):
                        b.a.charge += float(btype.q2)
                        b.b.charge += float(btype.q1)