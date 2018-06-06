# ******************************************************************************
# pysimm.forcefield.dreiding module
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


class Dreiding(Forcefield):
    """pysimm.forcefield.Dreiding

    Forcefield object with typing rules for Dreiding model.
    By default reads data file in forcefields subdirectory.

    Attributes:
        ff_name: dreiding
        pair_style: lj
        ff_class: 1
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   os.pardir, os.pardir, 'dat', 'forcefields', 'dreiding.json')
        Forcefield.__init__(self, db_file)
        self.name = 'dreiding'
        self.pair_style = 'lj'
        self.bond_style = 'harmonic'
        self.angle_style = 'harmonic'
        self.dihedral_style = 'harmonic'
        self.improper_style = 'umbrella'
        self.ff_class = '1'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Dreiding.assign_ptypes

        Dreiding specific particle typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds defined.
        *** use System.add_particle_bonding() to ensure this ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        s.pair_style = self.pair_style
        all_types = set()
        s.add_particle_bonding()
        for p in s.particles:
            p.bond_orders = [x.order for x in p.bonds]
            if len(set(p.bond_orders)) == 1 and p.bond_orders[0] is None:
                error_print('error: bond orders are not set')
            p.bond_elements = [x.a.elem if p is x.b else x.b.elem for x in
                               p.bonds]
            p.nbonds = len(p.bond_elements)
            if p.linker:
                p.nbonds += 1
        for p in s.particles:
            if p.elem == 'H':
                p.type_name = 'H_'
            elif p.elem == 'C':
                if p.bond_orders and (4 in p.bond_orders or 'A' in p.bond_orders):
                    p.type_name = 'C_R'
                elif p.nbonds == 4:
                    p.type_name = 'C_3'
                elif p.nbonds == 3:
                    p.type_name = 'C_2'
                elif p.nbonds == 2:
                    p.type_name = 'C_1'
                else:
                    print 'cant type particle %s' % p.tag
                    return p
            elif p.elem == 'N':
                if 4 in p.bond_orders or 'A' in p.bond_orders:
                    p.type_name = 'N_R'
                elif 2 in p.bond_orders:
                    p.type_name = 'N_2'
                elif 3 in p.bond_orders:
                    p.type_name = 'N_1'
                elif 1 in p.bond_orders:
                    for pb in p.bonded_to:
                        if pb.elem == 'C' and pb.nbonds == 3:
                            p.type_name = 'N_2'
                    if not p.type_name:
                        p.type_name = 'N_3'
                else:
                    print 'cant type particle %s' % p.tag
                    return p
            elif p.elem == 'O':
                if 4 in p.bond_orders or 'A' in p.bond_orders:
                    p.type_name = 'O_R'
                elif 2 in p.bond_orders:
                    p.type_name = 'O_2'
                elif 3 in p.bond_orders:
                    p.type_name = 'O_1'
                elif 1 in p.bond_orders and len(set(p.bond_orders)) == 1:
                    p.type_name = 'O_3'
                else:
                    print 'cant type particle %s' % p.tag
                    return p
            elif p.elem == 'F':
                p.type_name = 'F_'
            elif p.elem == 'P':
                p.type_name = 'P_3'
            elif p.elem == 'S':
                p.type_name = 'S_3'
            elif p.elem == 'Cl':
                p.type_name = 'Cl'
            else:
                print 'cant type particle %s' % p.tag
                return p
            all_types.add(self.particle_types.get(p.type_name)[0])

        for pt in all_types:
            pt = pt.copy()
            s.particle_types.add(pt)

        for p in s.particles:
            pt = s.particle_types.get(p.type_name)
            if pt:
                p.type = pt[0]

    def assign_btypes(self, s):
        """pysimm.forcefield.Dreiding.assign_btypes

        Dreiding specific bond typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds, type and type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        all_types = set()
        s.bond_style = self.bond_style
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
        """pysimm.forcefield.Dreiding.assign_atypes

        Dreiding specific angle typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds, type and type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        all_types = set()
        s.angle_style = self.angle_style
        s.add_particle_bonding()
        for p in s.particles:
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
        """pysimm.forcefield.Dreiding.assign_dtypes

        Dreiding specific dihedral typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds, type and type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        all_types = set()
        s.dihedral_style = self.dihedral_style
        for b in s.bonds:
            for p1 in b.a.bonded_to:
                for p2 in b.b.bonded_to:
                    if p1 is b.b or p2 is b.a:
                        continue
                    unique = True
                    for d in s.dihedrals:
                        if ((d.a == p1 and d.b == b.a and
                                d.c == b.b and d.d == p2) or
                            (d.a == p2 and d.b == b.b and
                                d.c == b.a and d.d == p1)):
                            unique = False
                    if unique:
                        if len(b.a.type.name) > 2 and b.a.type.name[2] == '1':
                            continue
                        if len(b.b.type.name) > 2 and b.b.type.name[2] == '1':
                            continue
                        p1_name = p1.type.name
                        a_name = b.a.type.name
                        b_name = b.b.type.name
                        p2_name = p2.type.name
                        dt = self.dihedral_types.get('%s,%s,%s,%s'
                                                     % (p1_name, a_name,
                                                        b_name, p2_name))
                        if dt:
                            if len(dt) == 1:
                                all_types.add(dt[0])
                                s.dihedrals.add(Dihedral(type_name=dt[0].name,
                                                         a=p1, b=b.a,
                                                         c=b.b, d=p2))
                            else:
                                index = 0
                                x = 5
                                for i in range(len(dt)):
                                    if dt[i].name.count('X') < x:
                                        index = i
                                        x = dt[i].name.count('X')
                                dt = dt[index]
                                all_types.add(dt)
                                s.dihedrals.add(Dihedral(type_name=dt.name,
                                                         a=p1, b=b.a,
                                                         c=b.b, d=p2))
                        else:
                            print ('I cant type this dihedral %s,%s,%s,%s'
                                   % (p1_name, a_name, b_name, p2_name))

        for dt in all_types:
            dt = dt.copy()
            s.dihedral_types.add(dt)

        for d in s.dihedrals:
            dt = s.dihedral_types.get(d.type_name)
            if dt:
                d.type = dt[0]

    def assign_itypes(self, s):
        """pysimm.forcefield.Dreiding.assign_itypes

        Dreiding specific improper typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds, type and type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        all_types = set()
        s.improper_style = self.improper_style
        for p in s.particles:
            if len(p.bonded_to) == 3:
                for perm in permutations(p.bonded_to, 3):
                    p1_name = perm[0].type.eq_improper or perm[0].type.name
                    p2_name = perm[1].type.eq_improper or perm[1].type.name
                    p3_name = perm[2].type.eq_improper or perm[2].type.name
                    it = self.improper_types.get(','.join([p.type.name, p1_name,
                                                           p2_name, p3_name]), order=True)
                    if it:
                        all_types.add(it[0])
                        bonded_to = p.bonded_to.get('all')
                        s.impropers.add(Improper(type_name=it[0].name,
                                                 a=p, b=bonded_to[0],
                                                 c=bonded_to[1],
                                                 d=bonded_to[2]))
                        break

        for it in all_types:
            it = it.copy()
            s.improper_types.add(it)

        for i in s.impropers:
            it = s.improper_types.get(i.type_name)
            if it:
                i.type = it[0]

    def assign_charges(self, s, charges='gasteiger'):
        """pysimm.forcefield.Dreiding.assign_charges

        Charge assignment. Gasteiger is default for now.

        Args:
            s: :class:`~pysimm.system.System`
            charges: gasteiger

        Returns:
            None
        """
        if charges == 'gasteiger':
            print('adding gasteiger charges')
            gasteiger.set_charges(s)