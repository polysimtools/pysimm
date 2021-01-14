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
import json
import os
import re
import sys
from itertools import permutations

import numpy

from . import gasteiger
from .. import error_print
from ..system import Angle, Dihedral, Improper, ParticleType
from .forcefield import Forcefield
from ..utils import ItemContainer


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
                os.pardir, 'data', 'forcefields', 'charmm.json'
            )
        Forcefield.__init__(self, db_file)

        with open(db_file) as f:
            j = json.loads(f.read())
        self.nondiag_lj_types = ItemContainer()
        for elem in j.get('nondiagonal_lj'):
            self.nondiag_lj_types.add(ParticleType(**elem))

        self.name = 'charmm'
        self.pair_style = 'lj/charmm'
        self.bond_style = 'harmonic'
        self.angle_style = 'charmm'
        self.dihedral_style = 'charmm'
        self.improper_style = 'harmonic'
        self.ff_class = '1'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Charmm.assign_ptypes

        Charmm specific particle typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds defined.
        *** use System.add_particle_bonding() to ensure this ***

        *** Not entirely inclusive - some atom types not used ***

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None
        """

        s.pair_style = self.pair_style
        s.add_particle_bonding()
        for p in s.particles:
            p.bond_orders = [x.order for x in p.bonds]
            if None in p.bond_orders:
                error_print('error: bond orders are not set')
            p.bond_elements = [x.a.elem if p is x.b else x.b.elem for x in p.bonds]
            p.nbonds = len(p.bond_elements)
            if p.linker:
                p.nbonds += 1

        for p in s.particles:
            if not p.type_name:

                if p.elem == 'C':
                    # General definition of a carbon in ethers
                    if (all(p.bond_orders) == 1) and (p.nbonds == 4):
                        rng_count = __detect_rings__(p, [5, 6])
                        if rng_count > 0 : # tetrahydrofuran (THF) or tetrahydropyran (THP)
                            p.type_name = 'CC32{}B'.format(rng_count)
                            for sb_p in p.bonded_to: # type all hydrogens connected to this atom
                                if sb_p.elem == 'H':
                                    sb_p.type_name = 'HCA25A'
                        else: # it is linear ether
                            hcount = len([True for l in p.bond_elements if l == 'H'])
                            p.type_name = 'CC3{}A'.format(hcount)
                            for sb_p in p.bonded_to: # type all hydrogens connected to this atom
                                if sb_p.elem == 'H':
                                    sb_p.type_name = 'HCA{}A'.format(hcount)

                if p.elem == 'O':
                    if (p.nbonds == 2) and (all(p.bond_orders) == 1) and all([t == 'C' for t in p.bond_elements]):
                        p.type_name = 'OC30A'
                        if __detect_rings__(p, [5, 6]):
                            p.type_name = 'OC305A'
        all_types = set()
        for p in s.particles:
            all_types.add(self.particle_types.get(p.type_name)[0])

        for pt in all_types:
            s.particle_types.add(pt.copy())

        for p in s.particles:
            pt = s.particle_types.get(p.type_name)
            if pt:
                p.type = pt[0]

        # Addition to normal FF setup: filling up the nondiagonal pair coefficient
        loc_lj_types = set()
        for p in s.particle_types:
            for p_ in s.particle_types:
                if p != p_:
                    atm_type = tuple(sorted([p.tag, p_.tag]))
                    if not(atm_type in [at.atm_types for at in loc_lj_types]):
                        tmp = self.nondiag_lj_types.get(','.join([p.name, p_.name]))
                        if len(tmp) > 0:
                            to_add = tmp[0].copy()
                            to_add.atm_types = atm_type
                            loc_lj_types.add(to_add)

        s.nondiag_lj_types = ItemContainer()
        for ljt in loc_lj_types:
            s.nondiag_lj_types.add(ljt)

    def assign_btypes(self, s):
        """pysimm.forcefield.Charmm.assign_btypes

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
                print('couldnt type this bond %s,%s'
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
        """pysimm.forcefield.Charmm.assign_atypes

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
        """pysimm.forcefield.Charmm.assign_dtypes

        CHARMM specific dihedral typing rules.
        Requires :class:`~pysimm.system.System` object :class:`~pysimm.system.Particle` objects have bonds, type
        and type.name defined.
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
                            print('I cant type this dihedral %s,%s,%s,%s'
                                  % (p1_name, a_name, b_name, p2_name))

        for dt in all_types:
            dt = dt.copy()
            dt.w = 0.0
            s.dihedral_types.add(dt)

        for d in s.dihedrals:
            dt = s.dihedral_types.get(d.type_name, item_wildcard=None)
            if dt:
                d.type = dt[0]

    def assign_itypes(self, s):
        """pysimm.forcefield.Charmm.assign_itypes

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
        """pysimm.forcefield.Charmm.assign_charges

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


def __parse_charmm__():
    import json

    kj2kcal = 4.184
    rounding = 8
    bnded_lib = 'ffbonded.itp'
    atmtype_lib = 'atomtypes.atp'
    nonbnded_lib = 'ffnonbonded.itp'

    DATA_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'data', 'forcefields', 'charmm')
    obj = {'angle_types': [], 'improper_types': [], 'bond_types': [], 'particle_types': [], 'dihedral_types': []}

    # Parsing bonded parameters of the FF
    try:
        with open(os.path.join(DATA_PATH, bnded_lib), 'r') as f:
            ff_file = f.readlines()
    except OSError:
        print('Required library file with CHARMM bonded parametrs \"{:}\" '
              'cannot be opened or read. \nExiting...'.format(bnded_lib))
        sys.exit(1)

    i = 0
    curr_type = ''
    while i < len(ff_file):
        line = ff_file[i].split()
        if ff_file[i][0] == '[':
            curr_type = line[1]
            if ff_file[i + 1].split()[1] == "'improper'":
                curr_type = 'impropertypes'
                i += 1
            print(curr_type)
            i += 1

        elif line:
            try:
                if curr_type == 'bondtypes':
                    k = round(float(line[4]) / (2 * kj2kcal * 100), rounding)
                    b = round(float(line[3]) * 10, rounding)
                    name = ','.join(line[0:2])
                    rname = ','.join(reversed(line[0:2]))
                    obj['bond_types'].append({'k': k, 'tag': name, 'r0': b, 'name': name, 'rname': rname})

                elif curr_type == 'angletypes':
                    theta0 = round(float(line[4]), rounding)
                    ktheta = round(float(line[5]) / (2 * kj2kcal), rounding)
                    ub0 = round(10 * float(line[6]), rounding)
                    kub = round(float(line[7]) / (2 * kj2kcal), rounding)
                    name = ','.join(line[0:3])
                    rname = ','.join(reversed(line[0:3]))
                    obj['angle_types'].append(
                        {'theta0': theta0, 'tag': name, 'k': ktheta, 'r_ub': ub0, 'k_ub': kub, 'name': name,
                         'rname': rname})

                elif curr_type == 'impropertypes':
                    k = round(float(line[6]) / (2 * kj2kcal), rounding)
                    x0 = round(float(line[5]), rounding)
                    name = ','.join(line[0:4])
                    rname = ','.join(reversed(line[0:4]))
                    obj['improper_types'].append({'k': k, 'tag': name, 'x0': x0, 'name': name, 'rname': rname})

                elif curr_type == 'dihedraltypes':
                    d = round(float(line[5]), rounding)
                    k = round(float(line[6]) / kj2kcal, rounding)
                    n = int(line[7])
                    name = ','.join(line[0:4])
                    rname = ','.join(reversed(line[0:4]))
                    obj['dihedral_types'].append({'d': d, 'k': k, 'tag': name, 'n': n, 'name': name, 'rname': rname})
            except ValueError:
                print('improper value at line', i)
            except IndexError:
                print('missing value at line', i)
        i += 1

    # Parsing non-bonded parameters of the FF
    with open('../data/elements_by_mass.json', 'r') as pntr:
        elemsDict = json.load(pntr)

    try:
        with open(os.path.join(DATA_PATH, nonbnded_lib), 'r') as nb_file:
            try:
                with open(os.path.join(DATA_PATH, atmtype_lib), 'r') as f:
                    atp_data = f.read()
            except OSError:
                print('Required library file with CHARMM atom types \"{:}\" '
                      'cannot be opened or read. \nExiting...'.format(atmtype_lib))
                sys.exit(1)

            fields = ['name', 'epsilon', 'sigma', 'elem', 'tag', 'mass', 'desc']
            for line in nb_file:
                if not (line[0] in [';', '#', '[', '\n']):
                    line = line.strip().split()
                    if len(line) > 6:

                        # checking validity of the description field
                        descr = re.findall('(?<=' + '{:>6}    '.format(line[0]) + '[\d| ][\d| ]\d\.\d{5} ; ).*', atp_data)
                        if len(descr) > 0:
                            descr = descr[0]
                        else:
                            descr = ''

                        # checking validity of the element name field
                        elemname = line[1]
                        if len(line[1]) > 0:
                            if int(line[1]) > 0:
                                elemname = elemsDict[line[1]]['symbol']

                        tmp = [line[0],
                               round(float(re.match('-?\d{1,}\.\d{1,}', line[6]).group(0)) / kj2kcal, rounding),
                               round(10 * float(re.match('-?\d{1,}\.\d{1,}', line[5]).group(0)), rounding),
                               elemname, line[0], float(line[2]), descr]
                        obj['particle_types'].append(dict(zip(fields, tmp)))

            # Parsing **non-diagonal** non-bonded parameters of the FF
            nb_file.seek(0)
            obj['nondiagonal_lj'] = []
            for line in nb_file:
                if not (line[0] in [';', '#', '[', '\n']):
                    line = line.strip().split()
                    if len(line) == 5:
                        i_name = line[0].strip()
                        j_name = line[1].strip()

                        obj['nondiagonal_lj'].append({'name': ','.join([i_name, j_name]),
                                                      'rname': ','.join([j_name, i_name]),
                                                      'epsilon': round(float(line[3]) / kj2kcal, rounding),
                                                      'sigma': round(10 * float(line[4]), rounding)
                                                      })
    except OSError:
        print('Required library file with CHARMM non-bonded parametrs \"{:}\" '
              'cannot be opened or read. \nExiting...'.format(nonbnded_lib))
        sys.exit(1)

    # Adding meta-information about FF styles and creating an output file
    chrm_type = Charmm()
    attr = ['pair_style', 'bond_style', 'angle_style', 'dihedral_style', 'improper_style']
    obj.update(dict(zip(attr, [getattr(chrm_type, t) for t in attr])))

    out_file = os.path.join(DATA_PATH, os.path.pardir, 'charmm.json')
    with open(out_file, 'w') as pntr:
        pntr.write(json.dumps(obj, indent=2))


def __detect_rings__(particle, orders):
    rng_count = 0
    neighb_list = []
    ordr_count = 2
    to_exclude = set([particle])

    neighb = []
    for p in to_exclude:
        neighb += [x.a if particle is x.b else x.b for x in p.bonds]

    while ordr_count < max(orders):
        to_include = []
        for p in neighb:
            to_include += [x.a if particle is x.b else x.b for x in p.bonds]

        neighb_list.append(set(to_include) - to_exclude)
        to_exclude = set(neighb)
        ordr_count += 1

    for o in orders:
        if particle in neighb_list[o - 3]:
            rng_count = o

    return rng_count