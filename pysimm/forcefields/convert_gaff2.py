#!/usr/bin/env python

##############################
#
#
#  forcefield.pcff Module
#
#
##############################

import sys
import math
from xml.dom import minidom
from xml.etree import ElementTree as ET
from pysimm.system import ParticleType, BondType, AngleType
from pysimm.system import DihedralType, ImproperType
from pysimm.forcefield import Gaff, element_names_by_mass
from pysimm.utils import Item, ItemContainer


def convert(file_, out):
    ff = Gaff(False)
    with file(file_) as f:
        for line in f:
            line = line.split()
            if len(line) > 0 and line[0] == 'ATOM_TYPES':
                print('reading atom types...')
                line = f.next().split()
                while line[0] != 'END':
                    desc = ' '.join(line[3:]).strip()
                    pt = ff.particle_types.get(line[0])
                    elem = element_names_by_mass.get(int(round(float(line[1]))))
                    if pt:
                        pt = pt[0]
                        pt.mass = float(line[1])
                        pt.elem = elem
                        pt.desc = desc
                    else:
                        ff.particle_types.add(ParticleType(mass=float(line[1]),
                                                           elem=elem,
                                                           name=line[0],
                                                           tag=line[0],
                                                           desc=desc))
                    line = f.next().split()

            if len(line) > 0 and line[0] == 'BOND_TYPES':
                print('reading bond types...')
                line = f.next()
                line = f.next()
                while line.split()[0] != 'END':
                    p1, p2 = map(lambda a: a.strip(), line[0:5].split('-'))
                    line = line[5:].split()
                    name = '%s,%s' % (p1, p2)
                    bt = ff.bond_types.get(name)
                    if bt:
                        bt = bt[0]
                        bt.k = float(line[0])
                        bt.r0 = float(line[1])
                    else:
                        ff.bond_types.add(BondType(tag=name,
                                                   name=name,
                                                   k=float(line[0]),
                                                   r0=float(line[1])))
                    line = f.next()

            if len(line) > 0 and line[0] == 'ANGLE_TYPES':
                print('reading angle types...')
                line = f.next()
                while line.split()[0] != 'END':
                    p1, p2, p3 = map(lambda a: a.strip(), line[0:8].split('-'))
                    name = ','.join([p1, p2, p3])
                    line = line[8:].split()
                    k = float(line[0])
                    theta0 = float(line[1])
                    bt = ff.angle_types.get(name)
                    if bt:
                        bt = bt[0]
                        bt.k = k
                        bt.theta0 = theta0
                    else:
                        ff.angle_types.add(AngleType(tag=name,
                                                     name=name,
                                                     k=k,
                                                     theta0=theta0))
                    line = f.next()

            if len(line) > 0 and line[0] == 'DIHEDRAL_TYPES':
                print('reading dihedral types...')
                line = f.next()
                while line.split()[0] != 'END':
                    if line[0] == '#':
                        line = f.next()
                        continue
                    p1, p2, p3, p4 = map(lambda a: a.strip(),
                                         line[0:11].split('-'))
                    name = ','.join([p1, p2, p3, p4])
                    line = line[11:].split()
                    k = float(line[1])/int(float(line[0]))
                    n = int(float(line[3]))
                    d = int(round(math.cos(float(line[2])/180*math.pi)))
                    dt = ff.dihedral_types.get(name, item_wildcard=None)
                    if dt:
                        dt = dt[0]
                        dt.m += 1
                        dt.k.append(k)
                        dt.n.append(n)
                        dt.d.append(d)
                    else:
                        ff.dihedral_types.add(DihedralType(tag=name,
                                                           name=name,
                                                           m=1,
                                                           k=[k],
                                                           n=[n],
                                                           d=[d]))
                    line = f.next()
                    
            if len(line) > 0 and line[0] == 'IMPROPER_TYPES':
                print('reading improper types...')
                line = f.next()
                while line.split()[0] != 'END':
                    if line[0] == '#':
                        line = f.next()
                        continue
                    p2, p3, p1, p4 = map(lambda a: a.strip(),
                                         line[0:11].split('-'))
                    name = ','.join([p1, p2, p3, p4])
                    line = line[11:].split()
                    k = float(line[0])
                    n = int(float(line[2]))
                    d = int(round(math.cos(float(line[1])/180*math.pi)))
                    dt = ff.dihedral_types.get(name)
                    if dt:
                        dt = dt[0]
                        dt.k = k
                        dt.n = n
                        dt.d = d
                    else:
                        ff.dihedral_types.add(DihedralType(tag=name,
                                                           name=name,
                                                           k=k,
                                                           n=n,
                                                           d=d))
                    line = f.next()

            if len(line) > 0 and line[0] == 'VDW_TYPES':
                print('reading nondonded terms...')
                line = f.next().split()
                while line[0] != 'END':
                    name = line[0]
                    sigma = float(line[1])*2/pow(2, 1/6.)
                    epsilon = float(line[2])
                    pt = ff.particle_types.get(name)
                    if pt:
                        pt = pt[0]
                        pt.sigma = sigma
                        pt.epsilon = epsilon
                    else:
                        ff.particle_types.add(ParticleType(tag=name,
                                                           name=name,
                                                           sigma=sigma,
                                                           epsilon=epsilon))
                    line = f.next().split()

    ff.write_json(out)


if __name__ == '__main__':
    convert(sys.argv[1], sys.argv[2])
