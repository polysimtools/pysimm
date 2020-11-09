# ******************************************************************************
# pysimm.forcefield.charmm module
# ******************************************************************************
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2020
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

from . import gasteiger
from ..system import Angle, Dihedral, Improper
from .forcefield import Forcefield


class Charmm(Forcefield):
    """pysimm.forcefield.Charmm

    Forcefield object with typing rules for CHARMM model.
    By default reads data file in forcefields subdirectory.

    Attributes:
        ff_name: charmm
        pair_style: lj
        ff_class: 1
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(
                os.path.dirname(
                    os.path.realpath(__file__)
                ),
                os.pardir, 'data', 'forcefields', 'gaff.json'
            )
        Forcefield.__init__(self, db_file)
        self.name = 'gaff'
        self.pair_style = 'lj'
        self.bond_style = 'harmonic'
        self.angle_style = 'harmonic'
        self.dihedral_style = 'fourier'
        self.improper_style = 'cvff'
        self.ff_class = '1'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Gaff.assign_ptypes

        Gaff specific particle typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds defined.
        *** use System.add_particle_bonding() to ensure this ***

        *** Not entirely inclusive - some atom types not used ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """
        pass

    def assign_btypes(self, s):
        """pysimm.forcefield.Gaff.assign_btypes

        Gaff specific bond typing rules.
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
        """pysimm.forcefield.Gaff.assign_atypes

        Gaff specific boanglend typing rules.
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
        """pysimm.forcefield.Gaff.assign_dtypes

        Gaff specific dihedral typing rules.
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
        """pysimm.forcefield.Gaff.assign_itypes

        Gaff specific improper typing rules.
        There are none.

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
        """pysimm.forcefield.Gaff.assign_charges

        Charge assignment. Gasteiger is default for now.

        Args:
            s: :class:`~pysimm.system.System`
            charges: gasteiger

        Returns:
            None
        """
        pass


def _parse_charmm():
    with open('/pysimm/pysimm/data/forcefields/charmm/ffbonded.itp') as f:
        ff_file = f.readlines()
    i = 0
    obj = {'angle_types':[], 'improper_types':[], 'bond_types':[], 'particle_types':[], 'dihedral_types':[]}
    curr_type = ''
    while i < len(ff_file):
        line = ff_file[i].split()
        if ff_file[i][0] == '[':
            curr_type = line[1]
            if ff_file[i+1].split()[1] == "'improper'":
                curr_type = 'impropertypes'
                i+=1
            print(curr_type)
            i+=1
            
        elif line != []:
            try:
                if curr_type == 'bondtypes':
                    k = float(line[4])/(2*4.18*100)
                    b = float(line[3])*10
                    name = ','.join(line[0:2])
                    rname = ','.join(reversed(line[0:2]))
                    obj['bond_types'].append({'k':k, 'tag':name, 'r0':b, 'name':name, 'rname':rname})

                elif curr_type == 'angletypes':
                    theta0 = float(line[4])
                    ktheta = float(line[5])/(2*4.18)
                    ub0 = float(line[6])
                    kub = float(line[7])
                    name = ','.join(line[0:3])
                    rname = ','.join(reversed(line[0:3]))
                    obj['angle_types'].append({'theta0':theta0, 'tag':name, 'k':ktheta, 'ub0': ub0, 'kub': kub, 'name':name, 'rname':rname})

                elif curr_type == 'impropertypes':
                    k = float(line[6])/(2*4.18)
                    x0 = float(line[5])
                    name = ','.join(line[0:4])
                    rname = ','.join(reversed(line[0:4]))
                    obj['improper_types'].append({'k':k, 'tag':name, 'x0':x0, 'name':name, 'rname':rname})
                
                elif curr_type == 'dihedraltypes':
                    d = float(line[5])
                    k = float(line[6])/(4.18)
                    n = int(line[7])
                    name = ','.join(line[0:4])
                    rname = ','.join(reversed(line[0:4]))
                    obj['dihedral_types'].append({'d':d, 'k':k, 'tag':name, 'n':n, 'name':name, 'rname':rname})
            except ValueError:
                print('improper value at line', i)
            except IndexError:
                print('missing value at line', i)
        i+=1