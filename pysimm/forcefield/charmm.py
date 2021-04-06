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
from .. import error_print, verbose_print
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
        self.dihedral_style = 'fourier'
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

        for b in s.bonds:
            if not b.order:
                b.order = 1

        for p in s.particles:
            p.bond_orders = [x.order for x in p.bonds]
            if None in p.bond_orders:
                error_print('error: bond orders are not set')
            p.bond_elements = [x.a.elem if p is x.b else x.b.elem for x in p.bonds]
            p.nbonds = len(p.bond_elements)
            if p.linker:
                p.nbonds += 1

        for p in s.particles:
            if p.elem == 'C':
                # Some general definition of an -sp3 carbons
                if (all(p.bond_orders) == 1) and (p.nbonds == 4):
                    n_partcls = [p_ for p_ in p.bonded_to if p_.elem == 'N']
                    if len(n_partcls) > 0 and (n_partcls[0].nbonds == 4) and (set(n_partcls[0].bond_elements) == {'C'}):
                        p.type_name = 'CG324'
                    else:
                        rng_count = __detect_rings__(p, [5, 6])
                        if rng_count == 0: # Linear sp3 carbon
                            hcount = p.bond_elements.count('H')
                            p.type_name = 'CG3{}1'.format(hcount)
                        elif rng_count > 0 : # tetrahydrofuran (THF) or tetrahydropyran (THP)
                            p.type_name = 'CG3C52'.format(rng_count)

                if 'A' in p.bond_orders:
                    p.type_name = 'CG2R61'

                if (p.nbonds == 3):  # carbonyl C condition
                    if set(p.bond_elements) == {'O', 'C', 'N'}:  # in amide
                        p.type_name = 'CG2O1'
                    if p.bond_elements.count('O') == 2:  # carbonyl C in esters or acids

                        tmp_part = [sb_p for sb_p in p.bonded_to if (sb_p.elem == 'O') and sb_p.nbonds == 2]
                        if len(tmp_part) > 0:  # deprotonated
                            p.type_name = 'CG2O2'
                        else:  # protonated
                            p.type_name = 'CG2O3'

                    if set(p.bond_elements) == {'O', 'C', 'H'}:  # carbonyl C in aldehyde
                        p.type_name = 'CG2O4'
                    if (p.bond_elements.count('O') == 1) and (p.bond_elements.count('C') == 2):  # in ketones
                        p.type_name = 'CC2O5'

            elif p.elem == 'O':
                if (p.nbonds == 2) and (all(p.bond_orders) == 1):  # ethers, esters
                    if p.bond_elements.count('C') == 2:
                        is_ester = False
                        for p_ in p.bonded_to:
                            if (p_.bond_elements.count('O') == 2) and (p_.nbonds == 3):
                                is_ester = True
                        if is_ester:
                            p.type_name = 'OG302'
                        else:
                            p.type_name = 'OG301'
                            rng_count = __detect_rings__(p, [5, 6])
                            if rng_count > 0:
                                p.type_name = 'OG3C{}1'.format(rng_count)

                if (p.nbonds == 1) and ('C' in p.bond_elements):  # sp2 oxygen
                    p_ = [t for t in p.bonded_to][0]
                    if set(p_.bond_elements) == {'O', 'C', 'N'}:  # in amide
                        p.type_name = 'OG2D1'
                    if p_.bond_elements.count('O') == 2:  # in acids
                        tmp_part = [sb_p for sb_p in p_.bonded_to if (sb_p.elem == 'O') and sb_p.nbonds == 2]
                        if len(tmp_part) > 0:
                            p.type_name = 'OG2D1'
                        else:
                            p.type_name = 'OG2D2'
                    if p_.bond_elements.count('C') == 2:  # in ketones
                        p.type_name = 'OG2D3'

                if ('S' in p.bond_elements) or ('P' in p.bond_elements):  # phosphate or sulfate
                    p.type_name = 'OG2P1'

                if (p.nbonds == 2) and (set(p.bond_elements) == {'C', 'H'}):
                    p_ = [t for t in p.bonded_to if t.elem != 'H'][0]
                    if p_.bond_elements.count('O') == 2:  # in acids
                        p.type_name = 'OG2D1'
                    if p_.bond_elements.count('O') == 1:  # hydroxyl oxygen
                        p.type_name = 'OG311'

                if(p.nbonds == 2) and all([t == 'H' for t in p.bond_elements]):  # water oxygen
                    p.type_name = 'OT'
                    for sb_p in p.bonded_to:  # type all hydrogens connected to this atom
                        if sb_p.elem == 'H':
                            sb_p.type_name = 'HT'

            elif p.elem == 'N':
                if (p.nbonds == 1) and ('C' in p.bond_elements):  # nitrile (or cyano) group
                    p.type_name = 'NG1T1'
                if (p.nbonds == 3) and (set(p.bond_elements) == {'H', 'N'}):  # hydrazine
                    p.type_name = 'NG3N1'
                if (p.nbonds == 3) and ('C' in p.bond_elements):  # amide
                    p.type_name = 'NG2S{}'.format(p.bond_elements.count('H'))
                if (p.nbonds == 4):
                    p.type_name = 'NG3P{}'.format(p.bond_elements.count('H'))

            elif p.elem == 'H':
                if p.bond_elements[0] == 'N':
                    p.type_name = 'HGP1'
                if p.bond_elements[0] == 'O':
                    p.type_name = 'HGP1'
                if p.bond_elements[0] == 'C':
                    host = [p_ for p_ in p.bonded_to][0]
                    nitrogen = [p_ for p_ in host.bonded_to if p_.elem == 'N']
                    if len(nitrogen) > 0 and (nitrogen[0].nbonds == 4):
                        p.type_name = 'HGP5'
                    else:
                        hcount = [pt for pt in p.bonded_to][0].bond_elements.count('H')
                        p.type_name = 'HGA{}'.format(hcount)

            elif p.elem == 'S':
                if p.nbonds == 4:
                    p.type_name = 'SG3O{}'.format(4 - p.bond_elements.count('O'))

            else:
                print('cant type particle %s' % p.tag)
                return p

        all_types = set()
        for p in s.particles:
            all_types.add(self.particle_types.get(p.type_name)[0])

        for pt in all_types:
            s.particle_types.add(pt.copy())

        for p in s.particles:
            pt = s.particle_types.get(p.type_name)
            if pt:
                p.type = pt[0]
        self.assign_extra_ljtypes(s)

    def assign_extra_ljtypes(self, s):
        """pysimm.forcefield.Charmm.assign_extra_ljtypes

                Addition to normal force field setup: filling up the non-diagonal interaction pair
                coefficients (coefficients for interaction of particles of different type).

                Assumes that all :class:`~pysimm.system.ParticleType` are defined for all particles in s

        Args:
            s: :class:`~pysimm.system.System`

        Returns:
            None

        """
        #
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

        if not s.nondiag_lj_types:
            s.nondiag_lj_types = ItemContainer()

        for ljt in loc_lj_types:
            if not s.nondiag_lj_types.get(ljt.name):
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
            dt.w = 1.0
            s.dihedral_types.add(dt)
        verbose_print('Dihedrals assigned successfully. \nIMPORTANT: all dihedral weighting factors '
                      '(coefficients to compensate for double counting in rings) are currently set to 1.0.\n'
                      'If those values are different for your system please multiply corresponding force constants '
                      'by the weights manually.\n')

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
            verbose_print('adding gasteiger charges')
            gasteiger.set_charges(s)


def __parse_charmm__():
    """
    Private method to read/convert CHARMM specific FF parameters definition files to .json file which PySIMM works with

    Returns:
        None

    """

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
        error_print('Required library file with CHARMM bonded parametrs \"{:}\" '
                    'cannot be opened or read. \nExiting...'.format(bnded_lib))
        sys.exit(1)

    i = 0
    curr_type = ''
    dig_names_check = []
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

                    if name not in dig_names_check:
                        obj['dihedral_types'].append({'tag': name, 'd': [d], 'k': [k], 'n': [n], 'm': 1,
                                                      'name': name, 'rname': rname})
                        dig_names_check.append(name)
                    else:
                        to_add = list(filter(lambda fields: fields['name'] == name, obj['dihedral_types']))[0]
                        to_add['d'].append(d)
                        to_add['k'].append(k)
                        to_add['n'].append(n)
                        to_add['m'] += 1

            except ValueError:
                print('improper value at line', i)
            except IndexError:
                print('missing value at line', i)
        i += 1

    # Parsing non-bonded parameters of the FF
    elems_json = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'data', 'elements_by_mass.json')
    with open(elems_json, 'r') as pntr:
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
    """
    Private method for analysing whether a given particle is a part of a ring structure

    Args:
        particle: :class:`~pysimm.system.Particle` reference
        orders: list of integers defining size of the rings which should be checked

    Returns:
        list of integers subset of orders which defines the sizes of the rings that contain particle;
        returns 0 if no cyclic structures of size orders are detected
    """
    rng_count = 0
    neighb_list = []
    ordr_count = 2
    to_exclude = {particle}

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
