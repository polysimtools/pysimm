# ******************************************************************************
# pysimm.forcefield.pcff module
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
from pysimm.utils import compare
from pysimm.system import Angle, Dihedral, Improper
from forcefield import Forcefield


class Pcff(Forcefield):
    """pysimm.forcefield.Pcff

    Forcefield object with typing rules for Pcff model.
    By default reads data file in forcefields subdirectory.

    Attributes:
        ff_name: pcff
        pair_style: class2
        ff_class: 2
        nb_mixing: sixth
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   os.pardir, os.pardir, 'dat', 'forcefields', 'pcff.xml')
        Forcefield.__init__(self, db_file)
        self.name = 'pcff'
        self.ff_class = '2'
        self.pair_style = 'class2'
        self.bond_style = 'class2'
        self.angle_style = 'class2'
        self.dihedral_style = 'class2'
        self.improper_style = 'class2'
        self.nb_mixing = 'sixth'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Pcff.assign_ptypes

        Pcff specific particle typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds defined.
        *** use System.add_particle_bonding() to ensure this ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        all_types = set()
        s.pair_style = self.pair_style
        s.add_particle_bonding()
        for p in s.particles:
            p.bond_elements = [x.a.elem if p is x.b else x.b.elem for x in
                               p.bonds]
            p.bond_orders = [x.order for x in p.bonds]
            p.nbonds = len(p.bond_elements)
            if p.linker:
                p.nbonds += 1
        for p in s.particles:
            if p.elem == 'H':
                if ('C' in p.bond_elements or 'Si' in p.bond_elements or
                        'H' in p.bond_elements or 'S' in p.bond_elements):
                    p.type_name = 'h'
                elif 'O' in p.bond_elements or 'N' in p.bond_elements:
                    p.type_name = 'h*'
                elif 'S' in p.bond_elements:
                    p.type_name = 'hs'
                else:
                    print 'dont think I can type this one'
                    return p
            elif p.elem == 'C':
                if p.nbonds == 4:
                    p.type_name = 'c'
                elif 3 in p.bond_orders:
                    p.type_name = 'ct'
                elif 4 in p.bond_orders or 'A' in p.bond_orders:
                    p.type_name = 'cp'
                elif p.nbonds == 3 and p.bond_elements.count('N') == 3:
                    p.type_name = 'cr'
                elif p.nbonds == 3 and p.bond_elements.count('O') == 3:
                    p.type_name = 'cz'
                elif (p.nbonds == 3 and p.bond_elements.count('O') == 2 and
                        'N' in p.bond_elements):
                    p.type_name = 'c_2'
                elif (p.nbonds == 3 and p.bond_elements.count('N') == 2 and
                        'O' in p.bond_elements):
                    p.type_name = 'c_2'
                elif (p.nbonds == 3 and p.bond_elements.count('O') == 2 and
                        'C' in p.bond_elements):
                    p.type_name = 'c_1'
                elif (p.nbonds == 3 and 'O' in p.bond_elements and
                        'C' in p.bond_elements and 'N' in p.bond_elements):
                    p.typ_name = 'c_1'
                elif (p.nbonds == 3 and p.bond_elements.count('C') == 2 and
                        'O' in p.bond_elements):
                    p.type_name = 'c_0'
                elif (p.nbonds == 3 and 'O' in p.bond_elements and
                        'C' in p.bond_elements and 'H' in p.bond_elements):
                    p.type_name = 'c_0'
                elif 2 in p.bond_orders:
                    p.type_name = 'c=2'
                else:
                    print 'dont think I can type this one'
                    return p
            elif p.elem == 'N':
                if 3 in p.bond_orders:
                    p.type_name = 'nt'
                elif 4 in p.bond_orders or 'A' in p.bond_orders:
                    p.type_name = 'nn'
                elif p.nbonds == 2 and p.bond_orders.count(1) == 2:
                    p.type_name = 'nz'
                elif p.nbonds == 3:
                    p.type_name = 'n'
                else:
                    print 'dont think I can type this one'
                    return p
            elif p.elem == 'O':
                ester = False
                if p.bond_elements.count('H') == 2:
                    p.type_name = 'o*'
                elif p.bond_elements.count('H') == 1:
                    p.type_name = 'o'
                elif p.nbonds == 2:
                    for bond in p.bonds:
                        temp_p = bond.a if p is bond.b else bond.b
                        if temp_p.nbonds != 4:
                            ester = True
                    if ester:
                        p.type_name = 'o_2'
                    else:
                        p.type_name = 'o'
                elif p.nbonds == 1:
                    bond = p.bonds.get('all')[0]
                    c = bond.a if p is bond.b else bond.b
                    if c.nbonds == 3 and c.elem == 'C':
                        p.type_name = 'o_1'
                    else:
                        p.type_name = 'o='
                else:
                    print 'dont think I can type this one'
                    return p
            elif p.elem == 'Cl':
                p.type_name = 'cl'
            elif p.elem == 'S':
                if p.nbonds == 1:
                    p.type_name == "s'"
                elif p.nbonds == 2:
                    if 'S' in p.bond_elements:
                        p.type_name = 's'
                    elif 'H' in p.bond_elements:
                        p.type_name = 's'
                    elif p.bond_elements.count('C') == 2:
                        p.type_name = 's'
                    else:
                        p.type_name = 's'
                elif p.nbonds == 4:
                    p.type_name = "s'"
            else:
                print 'dont think I can type this one'
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
        """pysimm.forcefield.Pcff.assign_btypes

        Pcff specific bond typing rules.
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
        """pysimm.forcefield.Pcff.assign_atypes

        Pcff specific angle typing rules.
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
                                if compare('%s,%s,%s' % (p1.type.name,
                                                         p.type.name,
                                                         p2.type.name),
                                           at[0].name, order=True):
                                    s.angles.add(Angle(type_name=at[0].name,
                                                       a=p1, b=p, c=p2))
                                    all_types.add(at[0])
                                elif compare('%s,%s,%s' % (p2.type.name,
                                                           p.type.name,
                                                           p1.type.name),
                                             at[0].name, order=True):
                                    s.angles.add(Angle(type_name=at[0].name,
                                                       a=p2, b=p, c=p1))
                                    all_types.add(at[0])
                                else:
                                    print('cannot distinguish between forward/'
                                          'backward for %s,%s,%s'
                                          % (p1.type.name,
                                             p.type.name,
                                             p2.type.name))
                            else:
                                print('I cant type this angle %s,%s,%s'
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
        """pysimm.forcefield.Pcff.assign_dtypes

        Pcff specific dihedral typing rules.
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
                        if ((d.a is p1 and d.b is b.a and
                             d.c is b.b and d.d is p2) or
                                (d.a is p2 and d.b is b.b and
                                 d.c is b.a and d.d is p1)):
                            unique = False
                    if unique:
                        p1_name = p1.type.eq_dihedral or p1.type.name
                        a_name = b.a.type.eq_dihedral or b.a.type.name
                        b_name = b.b.type.eq_dihedral or b.b.type.name
                        p2_name = p2.type.eq_dihedral or p2.type.name
                        dt = self.dihedral_types.get('%s,%s,%s,%s'
                                                     % (p1_name, a_name,
                                                        b_name, p2_name))
                        if dt:
                            if len(dt) == 1:
                                dt = dt[0]
                            else:
                                index = 0
                                x = 5
                                for i in range(len(dt)):
                                    if dt[i].name.count('X') < x:
                                        index = i
                                        x = dt[i].name.count('X')
                                dt = dt[index]
                            if compare('%s,%s,%s,%s' % (p1_name, a_name,
                                                        b_name, p2_name),
                                       dt.name, order=True):
                                all_types.add(dt)
                                s.dihedrals.add(Dihedral(type_name=dt.name,
                                                         a=p1, b=b.a,
                                                         c=b.b, d=p2))
                            elif compare('%s,%s,%s,%s' % (p2_name, b_name,
                                                          a_name, p1_name),
                                         dt.name, order=True):
                                all_types.add(dt)
                                s.dihedrals.add(Dihedral(type_name=dt.name,
                                                         a=p2, b=b.b,
                                                         c=b.a, d=p1))
                            else:
                                print('cannot distinguish between forward/'
                                      'backward for %s,%s,%s,%s'
                                      % (p1.type.name,
                                         b.a.type.name,
                                         b.b.type.name,
                                         p2.type.name))
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
        """pysimm.forcefield.Pcff.assign_itypes

        Pcff specific improper typing rules.
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
            for perm in permutations(p.bonded_to, 3):
                p1_name = perm[0].type.eq_improper or perm[0].type.name
                p2_name = perm[1].type.eq_improper or perm[1].type.name
                p3_name = perm[2].type.eq_improper or perm[2].type.name
                it = self.improper_types.get(','.join([p.type.name, p1_name,
                                                       p2_name, p3_name]), order=True)
                if it:
                    unique = True
                    for i in s.impropers:
                        if i.a is not p:
                            continue
                        if set([i.b, i.c, i.d]) == set([perm[0], perm[1],
                                                        perm[2]]):
                            unique = False
                            break
                    if unique:
                        all_types.add(it[0])
                        s.impropers.add(Improper(type_name=it[0].name,
                                                 a=p, b=perm[0], c=perm[1],
                                                 d=perm[2]))

        for it in all_types:
            it = it.copy()
            s.improper_types.add(it)

        for i in s.impropers:
            it = s.improper_types.get(i.type_name)
            if it:
                i.type = it[0]

    def assign_charges(self, s, charges='default'):
        """pysimm.forcefield.Pcff.assign_charges

        Default Pcff charge assignment. Gasteiger is also an option.

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
            print('adding default PCFF charges')
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