# ******************************************************************************
# pysimm.forcefield.cgenff module
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


class Cgenff(Forcefield):
    """pysimm.forcefield.Cgenff

    Forcefield object with typing rules for Cgenff model.
    By default reads data file in forcefields subdirectory.

    Attributes:
        ff_name: cgenff
        pair_style: lj
        bond_style: harmonic
        angle_style: harmonic
        dihedral_style: fourier
        improper_style: cvff
        ff_class: 1
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   os.pardir, os.pardir, 'dat', 'forcefields', 'cgenff.json')
        Forcefield.__init__(self, db_file)
        self.ff_name = 'cgenff'
        self.pair_style = 'lj'
        self.bond_style = 'harmonic'
        self.angle_style = 'harmonic'
        self.dihedral_style = 'fourier'
        self.improper_style = 'harmonic'
        self.ff_class = '1'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Gaff2.assign_ptypes

        Gaff2 specific particle typing rules.
        Requires System object Particle objects have Particle.bonds defined.
        *** use System.add_particle_bonding() to ensure this ***

        *** Not entirely inclusive - some atom types not used ***

        Args:
            s: pysimm.system.System

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
            if p.linker:
                p.nbonds += 1
        for p in s.particles:
            if p.elem == 'H':
                if 'O' in p.bond_elements:
                    p.type_name = 'HGP1'
                elif 'N' in p.bond_elements:
                    p.type_name = 'HGP1'
                elif 'S' in p.bond_elements:
                    p.type_name = 'HGP3'
                elif 'C' in p.bond_elements:
                    for pb in p.bonded_to:
                        if pb.elem == 'C':
                            if 2 in pb.bond_orders:
                                if pb.bond_elements.count('H') == 1:
                                    p.type_name = 'HGA4'
                                elif pb.bond_elements.count('H') == 2:
                                    p.type_name = 'HGA5'
                            elif 4 in pb.bond_orders or 'A' in pb.bond_orders:
                                p.type_name = 'HGR61'
                            elif pb.nbonds == 4:
                                if pb.bond_elements.count('F') == 1:
                                    p.type_name = 'HGA6'
                                elif pb.bond_elements.count('F') == 2:
                                    p.type_name = 'HGA7'
                                elif pb.bond_elements.count('H') == 1:
                                    p.type_name = 'HGA1'
                                elif pb.bond_elements.count('H') == 2:
                                    p.type_name = 'HGA2'
                                elif pb.bond_elements.count('H') == 3:
                                    p.type_name = 'HGA3'
                                elif pb.bond_elements.count('H') == 4:
                                    p.type_name = 'HGA3'
            elif p.elem == 'C':
                if p.nbonds == 3 and 'O' in p.bond_elements:
                    if p.bond_elements.count('N') == 1:
                        p.type_name = 'CG2O1'
                    elif p.bond_elements.count('O') == 1:
                        if p.bond_elements.count('H') >= 1:
                            p.type_name = 'CG2O4'
                        else:
                            p.type_name = 'CG2O5'
                    elif p.bond_elements.count('O') == 2:
                        p.type_name = 'CG2O2'
                    elif p.bond_elements.count('O') == 3 or p.bond_elements.count('N') == 2:
                        p.type_name = 'CG2O6'
                elif 4 in p.bond_orders or 'A' in p.bond_orders:
                    p.type_name = 'CG2R61'
                elif p.nbonds == 4:
                    if p.bond_elements.count('H') == 0:
                        if p.bond_elements.count('F') == 3:
                            p.type_name = 'CG302'
                        else:
                            p.type_name = 'CG301'
                    elif p.bond_elements.count('H') == 1:
                        if p.bond_elements.count('F') == 2:
                            p.type_name = 'CG312'
                        else:
                            p.type_name = 'CG311'
                    elif p.bond_elements.count('H') == 2:
                        if p.bond_elements.count('F') == 1:
                            p.type_name = 'CG322'
                        else:
                            p.type_name = 'CG321'
                    elif p.bond_elements.count('H') >= 3:
                        p.type_name = 'CG331'
                elif p.nbonds == 3 and not 'O' in p.bond_elements:
                        if p.bond_elements.count('H') == 2:
                            p.type_name = 'CG2DC3'
                        else:
                            p.type_name = 'CG2DC1'
                elif p.nbonds == 2:
                    if p.bond_elements.count('N') == 1:
                        p.type_name = 'CG1N1'
                    elif p.bond_elements.count('H') == 0:
                        p.type_name = 'CG1T1'
                    else:
                        p.type_name = 'CG1T2'
            elif p.elem == 'N':
                if 3 in p.bond_orders and p.bond_elements.count('C') == 1:
                    p.type_name = 'NG1T1'
                elif 2 in p.bond_orders and p.nbonds == 1 and p.bond_elements.count('N') == 1 and p.bonded_to[0].bond_elements.count('N') == 1:
                    p.type_name = 'NG1D1'
                elif p.nbonds <= 3:
                    amide = False
                    aromatic_ring = False
                    for pb in p.bonded_to:
                        if pb.elem == 'C':
                            if 4 in pb.bond_orders or 'A' in pb.bond_orders:
                                aromatic_ring = True
                            for b in pb.bonds:
                                bp = b.a if pb is b.b else b.b
                                if bp.elem == 'O' and b.order == 2:
                                    amide = True
                    if amide:
                        if p.bond_elements.count('H') == 2:
                            p.type_name = 'NG2S2'
                        else:
                            p.type_name = 'NG2S1'
                    elif aromatic_ring:
                        p.type_name = 'NG2R60'
                    elif nbonds == 3:
                        if p.bond_elements.count('C') == 3:
                            p.type_name = 'NG301'
                        elif p.bond_elements.count('C') == 2:
                            p.type_name = 'NG311'
                        elif p.bond_elements.count('C') == 1:
                            p.type_name = 'NG321'
                        elif p.bond_elements.count('H') == 3:
                            p.type_name == 'NG331'
                else:
                    print(p.elem, p.nbonds, p.bond_elements, p.bond_orders)
            elif p.elem == 'O':
                if p.nbonds == 1 and p.bond_elements.count('C') == 1:
                    pb = p.bonded_to[0]
                    if pb.bond_elements.count('H') == 1:
                        p.type_name = 'OG2D3'
                    else:
                        p.type_name = 'OG2D1'
                elif p.bond_elements.count('H') == 1:
                    p.type_name = 'OG311'
                elif p.nbonds == 2:
                    for pb in p.bonded_to:
                        if pb.elem == 'C':
                            if pb.bond_elements.count('O') == 2 and pb.bond_orders.count(2) == 1:
                                for b in pb.bonds:
                                    if b.a.elem == 'O' and b.b.elem == 'C' and b.order == 2:
                                        p.type_name = 'OG302'
                                    elif b.a.elem == 'C' and b.a.elem == 'O' and b.order == 2:
                                        p.type_name = 'OG302'
                    if not p.type_name:
                        p.type_name = 'OG301'
            elif p.elem == 'S':
                if p.nbonds == 1:
                    if p.bond_elements.count('C') == 2:
                        p.type_name = 'SG3O3'
                    else:
                        p.type_name = 'SG2D1'
                elif len(set(p.bond_orders)) == 1 and p.bond_orders[0] == 1 and p.nbonds == 2:
                    p.type_name = 'SG311'
            elif p.elem == 'Cl':
                pb = p.bonded_to[0]
                if pb.elem == 'C' and pb.bond_elements.count('Cl') == 1 or pb.bond_elements.count('Cl') == 2:
                    p.type_name = 'CLGA1'
                elif pb.elem == 'C' and pb.bond_elements.count('Cl') == 3:
                    p.type_name = 'CLGA3'
                elif pb.elem == 'C' and pb.bond_orders.count('A') == 1 or pb.bond_orders.count(4) == 1:
                    p.type_name = 'CLGR1'
            elif p.elem == 'Br':
                pb = p.bonded_to[0]
                if pb.elem == 'C' and pb.bond_elements.count('Br') == 1:
                    p.type_name == 'BRGA1'
                elif pb.elem == 'C' and pb.bond_elements.count('Br') == 2:
                    p.type_name == 'BRGA2'
                elif pb.elem == 'C' and pb.bond_elements.count('Br') == 3:
                    p.type_name == 'BRGA3'
                elif pb.elem == 'C' and pb.bond_orders.count('A') == 1 or pb.bond_orders.count(4) == 1:
                    p.type_name = 'BRGR1'
            elif p.elem == 'I':
                pb = p.bonded_to[0]
                if pb.elem == 'C' and pb.bond_orders.count('A') == 1 or pb.bond_orders.count(4) == 1:
                    p.type_name = 'IGR1'
            elif p.elem == 'F':
                pb = p.bonded_to[0]
                if pb.elem == 'C' and pb.bond_elements.count('F') == 1:
                    p.type_name == 'FGA1'
                elif pb.elem == 'C' and pb.bond_elements.count('F') == 2:
                    p.type_name == 'FGA2'
                elif pb.elem == 'C' and pb.bond_elements.count('F') == 3:
                    p.type_name == 'FGA3'
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
        """pysimm.forcefield.Gaff2.assign_btypes

        Gaff2 specific bond typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        all_types = set()
        s.bond_style = self.bond_style
        for b in s.bonds:
            if b.a.type.name == 'CG2DC1' and b.b.type.name == 'CG2DC1' and b.order == 1:
                bt = self.bond_types.get('%s,%s' % (b.a.type.name, 'CG2DC2'))
            else:
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
        """pysimm.forcefield.Gaff2.assign_atypes

        Gaff2 specific boanglend typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        all_types = set()
        s.angle_style = self.angle_style
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
        """pysimm.forcefield.Gaff2.assign_dtypes

        Gaff2 specific dihedral typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
            dt = s.dihedral_types.get(d.type_name, item_wildcard=None)
            if dt:
                d.type = dt[0]

    def assign_itypes(self, s):
        """pysimm.forcefield.Gaff2.assign_itypes

        Gaff2 specific improper typing rules.
        There are none.

        Args:
            s: pysimm.system.System

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
                        s.impropers.add(Improper(type_name=it[0].name,
                                                 a=p, b=p.bonded_to[0],
                                                 c=p.bonded_to[1],
                                                 d=p.bonded_to[2]))
                        break

        for it in all_types:
            it = it.copy()
            s.improper_types.add(it)

        for i in s.impropers:
            it = s.improper_types.get(i.type_name)
            if it:
                i.type = it[0]

    def assign_charges(self, s, charges='gasteiger'):
        """pysimm.forcefield.Gaff.assign_charges

        Charge assignment. Gasteiger is default for now.

        Args:
            s: pysimm.system.System
            charges: gasteiger

        Returns:
            None
        """
        if charges == 'gasteiger':
            print('adding gasteiger charges')
            gasteiger.set_charges(s)