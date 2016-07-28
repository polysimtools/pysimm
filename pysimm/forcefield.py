# ******************************************************************************
# pysimm.forcefield module
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

import os
import json
from xml.dom import minidom
from xml.etree import ElementTree as Et
from itertools import permutations, combinations

from pysimm import error_print
from pysimm import warning_print
from pysimm import debug_print
from pysimm import gasteiger
from pysimm.utils import PysimmError, Item, ItemContainer, compare
from pysimm.system import ParticleType, BondType, AngleType
from pysimm.system import Angle, Dihedral, Improper
from pysimm.system import DihedralType, ImproperType


element_names_by_mass = {1: 'H', 4: 'He', 7: 'Li', 9: 'Be', 11: 'B', 12: 'C',
                         14: 'N', 16: 'O', 19: 'F', 20: 'Ne', 23: 'Na',
                         24: 'Mg', 27: 'Al', 28: 'Si', 31: 'P', 32: 'S',
                         35: 'Cl', 39: 'K', 40: 'Ca', 80: 'Br', 127: 'I'}


class Forcefield(object):
    """pysimm.forcefield.Forcefield

    Base Forcefield class definition.
    Initialize with force field xml file.

    Attributes:
        ff_class: force field class (1 or 2)
        ff_name: force field name
        particle_types: pysimm.utils.ItemContainer for particle_types
        bond_types: pysimm.utils.ItemContainer for bond_types
        angle_types: pysimm.utils.ItemContainer for angle_types
        dihedral_types: pysimm.utils.ItemContainer for dihedral_types
        improper_types: pysimm.utils.ItemContainer for improper_types
    """
    def __init__(self, file_=None):
        self.ff_class = 0
        self.ff_name = ''
        self.pair_style = None
        self.particle_types = ItemContainer()
        self.bond_style = None
        self.bond_types = ItemContainer()
        self.angle_style = None
        self.angle_types = ItemContainer()
        self.dihedral_style = None
        self.dihedral_types = ItemContainer()
        self.improper_style = None
        self.improper_types = ItemContainer()

        if not file_:
            return

        tree = Et.parse(file_)
        root = tree.getroot()
        ptypes = root.find('ParticleTypes')
        btypes = root.find('BondTypes')
        atypes = root.find('AngleTypes')
        dtypes = root.find('DihedralTypes')
        itypes = root.find('ImproperTypes')

        for ptype in ptypes:
            pt = ParticleType()
            for k, v in ptype.attrib.items():
                try:
                    v = float(v)
                except ValueError:
                    pass
                setattr(pt, k, v)
            pt.tag = pt.name
            self.particle_types.add(pt)

        for btype in btypes:
            bt = BondType()
            for k, v in btype.attrib.items():
                try:
                    v = float(v)
                except ValueError:
                    pass
                setattr(bt, k, v)
            bt.tag = bt.name
            bt.rname = ','.join(reversed(bt.name.split(',')))
            self.bond_types.add(bt)

        for atype in atypes:
            at = AngleType()
            for k, v in atype.attrib.items():
                try:
                    v = float(v)
                except ValueError:
                    pass
                setattr(at, k, v)
            at.tag = at.name
            at.rname = ','.join(reversed(at.name.split(',')))
            self.angle_types.add(at)

        for dtype in dtypes:
            dt = DihedralType()
            for k, v in dtype.attrib.items():
                if k != 'd' and k != 'n':
                    try:
                        v = float(v)
                    except ValueError:
                        pass
                elif k == 'd' or k == 'n':
                    try:
                        v = int(v)
                    except ValueError:
                        pass
                setattr(dt, k, v)
            dt.tag = dt.name
            dt.rname = ','.join(reversed(dt.name.split(',')))
            self.dihedral_types.add(dt)

        for itype in itypes:
            it = ImproperType()
            for k, v in itype.attrib.items():
                try:
                    v = float(v)
                except ValueError:
                    pass
                setattr(it, k, v)
            it.tag = it.name
            it.rname = ','.join(reversed(it.name.split(',')))
            self.improper_types.add(it)
            
    
    def from_json(self, json_file):
        with file(json_file) as f:
            j = json.loads(f.read())
        self.ff_name = j.get('ff_name')
        self.ff_class = j.get('ff_class')
        self.pair_style = j.get('pair_style')
        self.bond_style = j.get('bond_style')
        self.angle_style = j.get('angle_style')
        self.dihedral_style = j.get('dihedral_style')
        self.improper_style = j.get('improper_style')
        
        for pt in j.get('particle_types'):
            self.particle_types.add(ParticleType(**pt))
        
        for bt in j.get('bond_types'):
            self.bond_types.add(BondType(**bt))
        
        for at in j.get('angle_types'):
            self.angle_types.add(AngleType(**at))
        
        for dt in j.get('dihedral_types'):
            self.dihedral_types.add(DihedralType(**dt))
        
        for it in j.get('improper_types'):
            self.improper_types.add(ImproperType(**it))
        
            
    def write_json(self, out):
        f = {}
        f['ff_name'] = self.ff_name
        f['ff_class'] = self.ff_class
        f['pair_style'] = self.pair_style
        f['bond_style'] = self.bond_style
        f['angle_style'] = self.angle_style
        f['dihedral_style'] = self.dihedral_style
        f['improper_style'] = self.improper_style
        f['particle_types'] = [vars(pt) for pt in self.particle_types]
        f['bond_types'] = [vars(bt) for bt in self.bond_types]
        f['angle_types'] = [vars(at) for at in self.angle_types]
        f['dihedral_types'] = [vars(dt) for dt in self.dihedral_types]
        f['improper_types'] = [vars(it) for it in self.improper_types]
        with file(out, 'w') as _file:
            _file.write(json.dumps(f, indent=4))
            

    def write(self, out):
        """pysimm.forcefield.Forcefield.write

        Write Forcefield object to xml format.

        Args:
            out: file name to write

        Returns:
            None
        """

        ff_elem = Et.Element('Forcefield')
        ff_elem.set('name', self.ff_name)
        ff_elem.set('class', self.ff_class)

        ptypes = Et.SubElement(ff_elem, 'ParticleTypes')
        for pt in self.particle_types:
            ptype = Et.SubElement(ptypes, 'Type')
            for k, v in vars(pt).items():
                ptype.set(k, str(v))

        btypes = Et.SubElement(ff_elem, 'BondTypes')
        for bt in self.bond_types:
            btype = Et.SubElement(btypes, 'Type')
            for k, v in vars(bt).items():
                btype.set(k, str(v))

        atypes = Et.SubElement(ff_elem, 'AngleTypes')
        for at in self.angle_types:
            atype = Et.SubElement(atypes, 'Type')
            for k, v in vars(at).items():
                atype.set(k, str(v))

        dtypes = Et.SubElement(ff_elem, 'DihedralTypes')
        for dt in self.dihedral_types:
            dtype = Et.SubElement(dtypes, 'Type')
            for k, v in vars(dt).items():
                dtype.set(k, str(v))

        itypes = Et.SubElement(ff_elem, 'ImproperTypes')
        for it in self.improper_types:
            itype = Et.SubElement(itypes, 'Type')
            for k, v in vars(it).items():
                itype.set(k, str(v))

        with open(out, 'w') as f:
            f.write(
                minidom.parseString(Et.tostring(ff_elem, 'utf-8')).toprettyxml(
                    indent="  "))


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
                                   'forcefields', 'tip3p.xml')
        Forcefield.__init__(self, db_file)
        self.ff_name = 'tip3p'
        self.pair_style = 'lj'
        self.ff_class = '1'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Tip3p.assign_ptypes

        Tip3p specific particle typing rules.
        Requires System object Particle objects have Particle.bonds defined.
        *** use System.add_particle_bonding() to ensure this ***

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
        Requires System object Particle objects have Particle.type and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
            s: pysimm.system.System

        Returns:
            None
        """
        pass

    def assign_itypes(self, s):
        """pysimm.forcefield.Tip3p.assign_itypes

        Tip3p specific improper typing rules.
        There are none.

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        pass

    def assign_charges(self, s, charges='default'):
        """pysimm.forcefield.Tip3p.assign_charges

        Tip3p specific charge assignment.
        There are none.

        Args:
            s: pysimm.system.System
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


class Trappe(Forcefield):
    """pysimm.forcefield.Trappe

    Forcefield object with conversion rules for Trappe model.
    By default reads data file in forcefields subdirectory.

    *** IN PROGRESS - CONVERSION RULES NOT DEFINED ***

    Attributes:
        ff_name: trappe
        pair_style: lj
        ff_class: 1
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'forcefields', 'trappe.xml')
        Forcefield.__init__(self, db_file)
        self.ff_name = 'trappe'
        self.pair_style = 'lj'
        self.ff_class = '1'


class Gaff(Forcefield):
    """pysimm.forcefield.Gaff

    Forcefield object with typing rules for Gaff model.
    By default reads data file in forcefields subdirectory.

    Attributes:
        ff_name: gaff
        pair_style: lj
        ff_class: 1
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'forcefields', 'gaff.xml')
        Forcefield.__init__(self, db_file)
        self.ff_name = 'gaff'
        self.pair_style = 'lj'
        self.ff_class = '1'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Gaff.assign_ptypes

        Gaff specific particle typing rules.
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
                    water = False
                    for pb in p.bonded_to:
                        if pb.elem == 'O' and pb.bond_elements.count('H') == 2:
                            water = True
                    if water:
                        p.type_name = 'hw'
                    else:
                        p.type_name = 'ho'
                elif 'N' in p.bond_elements:
                    p.type_name = 'hn'
                elif 'P' in p.bond_elements:
                    p.type_name = 'hp'
                elif 'S' in p.bond_elements:
                    p.type_name = 'hs'
                elif 'C' in p.bond_elements:
                    for pb in p.bonded_to:
                        if pb.elem == 'C':
                            elctrwd = (pb.bond_elements.count('N') +
                                       pb.bond_elements.count('O') +
                                       pb.bond_elements.count('F') +
                                       pb.bond_elements.count('P') +
                                       pb.bond_elements.count('S') +
                                       pb.bond_elements.count('Cl'))
                            if 4 in pb.bond_orders or 'A' in pb.bond_orders:
                                p.type_name = 'ha'
                            elif elctrwd == 0:
                                p.type_name = 'hc'
                            elif pb.nbonds == 4 and elctrwd == 1:
                                p.type_name = 'h1'
                            elif pb.nbonds == 4 and elctrwd == 2:
                                p.type_name = 'h2'
                            elif pb.nbonds == 4 and elctrwd == 3:
                                p.type_name = 'h3'
                            elif pb.nbonds == 3 and elctrwd == 1:
                                p.type_name = 'h4'
                            elif pb.nbonds == 3 and elctrwd == 2:
                                p.type_name = 'h5'
            elif p.elem == 'C':
                if p.nbonds == 3 and 'O' in p.bond_elements:
                    p.type_name = 'c'
                elif 4 in p.bond_orders or 'A' in p.bond_orders:
                    for pb in p.bonded_to:
                        if pb.elem != 'C' and pb.elem != 'H':
                            p.type_name = 'cc'
                    if not p.type_name:
                        p.type_name = 'ca'
                elif p.nbonds == 4:
                    p.type_name = 'c3'
                elif p.nbonds == 3 and not 'O' in p.bond_elements:
                        p.type_name = 'c2'
                elif p.nbonds == 2:
                    p.type_name = 'c1'
            elif p.elem == 'F':
                p.type_name = 'f'
            elif p.elem == 'Cl':
                p.type_name = 'cl'
            elif p.elem == 'Br':
                p.type_name = 'br'
            elif p.elem == 'I':
                p.type_name = 'i'
            elif p.elem == 'N':
                if 2 in p.bond_orders and p.nbonds == 2:
                    p.type_name = 'n2'
                elif 2 in p.bond_orders and p.nbonds == 3:
                    p.type_name = 'na'
                elif 3 in p.bond_orders:
                    p.type_name = 'n1'
                elif p.bond_elements.count('O') == 2:
                    p.type_name = 'no'
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
                        p.type_name = 'n'
                    elif (4 in p.bond_orders or 'A' in p.bond_orders) and 'C' in p.bond_elements:
                        p.type_name = 'nc'
                    elif aromatic_ring:
                        p.type_name = 'nh'
                    else:
                        p.type_name = 'n3'
                elif p.nbonds == 4:
                    p.type_name = 'n4'
                else:
                    print(p.elem, p.nbonds, p.bond_elements, p.bond_orders)
            elif p.elem == 'O':
                if p.nbonds == 1:
                    p.type_name = 'o'
                elif p.bond_elements.count('H') == 2:
                    p.type_name = 'ow'
                elif p.bond_elements.count('H') == 1:
                    p.type_name = 'oh'
                else:
                    p.type_name = 'os'
            elif p.elem == 'P':
                if 4 in p.bond_orders or 'A' in p.bond_orders:
                    p.type_name = 'pc'
                elif p.nbonds == 2:
                    p.type_name = 'p2'
                elif p.nbonds == 3 and p.bond_elements.count('H') == 3:
                    p.type_name = 'p3'
                elif p.nbonds == 3:
                    p.type_name = 'p4'
                elif p.nbonds == 4:
                    p.type_name = 'p5'
            elif p.elem == 'S':
                if p.nbonds == 1:
                    p.type_name = 's'
                elif p.nbonds == 2 and 2 in p.bond_orders:
                    p.type_name = 's2'
                elif p.nbonds == 3:
                    p.type_name = 's4'
                elif p.nbonds == 4:
                    p.type_name = 's6'
                elif len(set(p.bond_orders)) == 1 and p.bond_orders[0] == 1:
                    if 'H' in p.bond_elements:
                        p.type_name = 'sh'
                    elif p.nbonds == 2:
                        p.type_name = 'ss'
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
        """pysimm.forcefield.Gaff.assign_btypes

        Gaff specific bond typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
        """pysimm.forcefield.Gaff.assign_atypes

        Gaff specific boanglend typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
        """pysimm.forcefield.Gaff.assign_dtypes

        Gaff specific dihedral typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        all_types = set()
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
                                                        b_name, p2_name),
                                                     item_wildcard=None)
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
            s: pysimm.system.System

        Returns:
            None
        """
        pass

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


class Gaff2(Forcefield):
    """pysimm.forcefield.Gaff

    Forcefield object with typing rules for Gaff model.
    By default reads data file in forcefields subdirectory.

    Attributes:
        ff_name: gaff
        pair_style: lj
        ff_class: 1
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'forcefields', 'gaff.xml')
        Forcefield.__init__(self, db_file)
        self.ff_name = 'gaff2'
        self.pair_style = 'lj'
        self.bond_style = ' harmonic'
        self.angle_style = 'harmonic'
        self.dihedral_style = 'fourier'
        self.improper_style = 'cvff'
        self.ff_class = '1'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Gaff.assign_ptypes

        Gaff specific particle typing rules.
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
                    water = False
                    for pb in p.bonded_to:
                        if pb.elem == 'O' and pb.bond_elements.count('H') == 2:
                            water = True
                    if water:
                        p.type_name = 'hw'
                    else:
                        p.type_name = 'ho'
                elif 'N' in p.bond_elements:
                    p.type_name = 'hn'
                elif 'P' in p.bond_elements:
                    p.type_name = 'hp'
                elif 'S' in p.bond_elements:
                    p.type_name = 'hs'
                elif 'C' in p.bond_elements:
                    for pb in p.bonded_to:
                        if pb.elem == 'C':
                            elctrwd = (pb.bond_elements.count('N') +
                                       pb.bond_elements.count('O') +
                                       pb.bond_elements.count('F') +
                                       pb.bond_elements.count('P') +
                                       pb.bond_elements.count('S') +
                                       pb.bond_elements.count('Cl'))
                            if 4 in pb.bond_orders or 'A' in pb.bond_orders:
                                p.type_name = 'ha'
                            elif elctrwd == 0:
                                p.type_name = 'hc'
                            elif pb.nbonds == 4 and elctrwd == 1:
                                p.type_name = 'h1'
                            elif pb.nbonds == 4 and elctrwd == 2:
                                p.type_name = 'h2'
                            elif pb.nbonds == 4 and elctrwd == 3:
                                p.type_name = 'h3'
                            elif pb.nbonds == 3 and elctrwd == 1:
                                p.type_name = 'h4'
                            elif pb.nbonds == 3 and elctrwd == 2:
                                p.type_name = 'h5'
            elif p.elem == 'C':
                if p.nbonds == 3 and 'O' in p.bond_elements:
                    p.type_name = 'c'
                elif 4 in p.bond_orders or 'A' in p.bond_orders:
                    for pb in p.bonded_to:
                        if pb.elem != 'C' and pb.elem != 'H':
                            p.type_name = 'cc'
                    if not p.type_name:
                        p.type_name = 'ca'
                elif p.nbonds == 4:
                    p.type_name = 'c3'
                elif p.nbonds == 3 and not 'O' in p.bond_elements:
                        p.type_name = 'c2'
                elif p.nbonds == 2:
                    p.type_name = 'c1'
            elif p.elem == 'F':
                p.type_name = 'f'
            elif p.elem == 'Cl':
                p.type_name = 'cl'
            elif p.elem == 'Br':
                p.type_name = 'br'
            elif p.elem == 'I':
                p.type_name = 'i'
            elif p.elem == 'N':
                if 2 in p.bond_orders and p.nbonds == 2:
                    p.type_name = 'n2'
                elif 2 in p.bond_orders and p.nbonds == 3:
                    p.type_name = 'na'
                elif 3 in p.bond_orders:
                    p.type_name = 'n1'
                elif p.bond_elements.count('O') == 2:
                    p.type_name = 'no'
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
                        p.type_name = 'n'
                    elif (4 in p.bond_orders or 'A' in p.bond_orders) and 'C' in p.bond_elements:
                        p.type_name = 'nc'
                    elif aromatic_ring:
                        p.type_name = 'nh'
                    else:
                        p.type_name = 'n3'
                elif p.nbonds == 4:
                    p.type_name = 'n4'
                else:
                    print(p.elem, p.nbonds, p.bond_elements, p.bond_orders)
            elif p.elem == 'O':
                if p.nbonds == 1:
                    p.type_name = 'o'
                elif p.bond_elements.count('H') == 2:
                    p.type_name = 'ow'
                elif p.bond_elements.count('H') == 1:
                    p.type_name = 'oh'
                else:
                    p.type_name = 'os'
            elif p.elem == 'P':
                if 4 in p.bond_orders or 'A' in p.bond_orders:
                    p.type_name = 'pc'
                elif p.nbonds == 2:
                    p.type_name = 'p2'
                elif p.nbonds == 3 and p.bond_elements.count('H') == 3:
                    p.type_name = 'p3'
                elif p.nbonds == 3:
                    p.type_name = 'p4'
                elif p.nbonds == 4:
                    p.type_name = 'p5'
            elif p.elem == 'S':
                if p.nbonds == 1:
                    p.type_name = 's'
                elif p.nbonds == 2 and 2 in p.bond_orders:
                    p.type_name = 's2'
                elif p.nbonds == 3:
                    p.type_name = 's4'
                elif p.nbonds == 4:
                    p.type_name = 's6'
                elif len(set(p.bond_orders)) == 1 and p.bond_orders[0] == 1:
                    if 'H' in p.bond_elements:
                        p.type_name = 'sh'
                    elif p.nbonds == 2:
                        p.type_name = 'ss'
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
        """pysimm.forcefield.Gaff.assign_btypes

        Gaff specific bond typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
        """pysimm.forcefield.Gaff.assign_atypes

        Gaff specific boanglend typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
        """pysimm.forcefield.Gaff.assign_dtypes

        Gaff specific dihedral typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        all_types = set()
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
                                                        b_name, p2_name),
                                                     item_wildcard=None)
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
            s: pysimm.system.System

        Returns:
            None
        """
        pass

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


class Cgenff(Forcefield):
    """pysimm.forcefield.Cgenff

    (Currently) empty Forcefield object for Cgenff model.
    By default reads data file in forcefields subdirectory.

    Attributes:
        ff_name: cgenff
        ff_class: 1
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'forcefields', 'cgenff.xml')
        Forcefield.__init__(self, db_file)
        self.ff_name = 'cgenff'
        self.ff_class = '1'


class Dreiding(Forcefield):
    """pysimm.forcefield.Dreiding

    Forcefield object with typing rules for Dreiding model.
    By default reads data file in forcefields subdirectory.

    Attributes:
        ff_name: dreiding
        pair_style: buck
        ff_class: 1
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'forcefields', 'dreiding.xml')
        Forcefield.__init__(self, db_file)
        self.ff_name = 'dreiding'
        self.pair_style = 'buck'
        self.ff_class = '1'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Dreiding.assign_ptypes

        Dreiding specific particle typing rules.
        Requires System object Particle objects have Particle.bonds defined.
        *** use System.add_particle_bonding() to ensure this ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        s.pair_style = self.pair_style
        all_types = set()
        for p in s.particles:
            p.bonded_to = [x.a if p is x.b else x.b for x in p.bonds]
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
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
        """pysimm.forcefield.Dreiding.assign_atypes

        Dreiding specific angle typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
        """pysimm.forcefield.Dreiding.assign_dtypes

        Dreiding specific dihedral typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        all_types = set()
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
                    if (unique and (b.a.type.name[2] != '1' and
                                    b.b.type.name[2] != '1')):
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
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        all_types = set()
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
        """pysimm.forcefield.Dreiding.assign_charges

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
                                   'forcefields', 'pcff.xml')
        Forcefield.__init__(self, db_file)
        self.ff_name = 'pcff'
        self.ff_class = '2'
        self.pair_style = 'class2'
        self.nb_mixing = 'sixth'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Pcff.assign_ptypes

        Pcff specific particle typing rules.
        Requires System object Particle objects have Particle.bonds defined.
        *** use System.add_particle_bonding() to ensure this ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        all_types = set()
        s.pair_style = self.pair_style
        for p in s.particles:
            p.bond_orders = [x.order for x in p.bonds]
            p.bond_elements = [x.a.elem if p is x.b else x.b.elem for x in
                               p.bonds]
            p.nbonds = len(p.bond_elements)
            if p.linker:
                p.nbonds += 1
        for p in s.particles:
            if p.elem == 'H':
                if ('C' in p.bond_elements or 'Si' in p.bond_elements or
                        'H' in p.bond_elements):
                    p.type_name = 'h'
                elif 'O' in p.bond_elements or 'N' in p.bond_elements:
                    p.type_name = 'h*'
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
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
        """pysimm.forcefield.Pcff.assign_atypes

        Pcff specific angle typing rules.
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

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
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        all_types = set()
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
        Requires System object Particle objects have Particle.bonds, Particle.type
        and Particle.type.name defined.
        *** use after assign_ptypes ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        all_types = set()
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
            s: pysimm.system.System
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

