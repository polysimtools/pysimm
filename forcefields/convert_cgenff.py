#!/usr/bin/env python

##############################
#
#
#  convert_cgenff script
#
#
##############################

import sys
from math import cos
from pysimm import system, forcefield


def convert(fname, out):
    ff = forcefield.Cgenff(False)
    with file(fname) as f:
        for line in f:
            if line.strip() == 'ATOMS':
                line = f.next()
                while line.strip() != 'BONDS':
                    if line.strip() and line[0] != '!':
                        data, desc = line.split('!')[0], '!'.join(line.split('!')[1:]).strip()
                        data = data.split()
                        ff.particle_types.add(system.ParticleType(name=data[2], tag=data[2],
                                                                  mass=float(data[3]),
                                                                  elem=data[2].split('G')[0].title(),
                                                                  desc=desc))
                    line = f.next()
            if line.strip() == 'BONDS':
                line = f.next()
                while line.strip() != 'ANGLES':
                    if line.strip() and line[0] != '!':
                        data, desc = line.split('!')[0], '!'.join(line.split('!')[1:]).strip()
                        data = data.split()
                        name = '{},{}'.format(data[0], data[1])
                        ff.bond_types.add(system.BondType(name=name, tag=name, k=float(data[2]),
                                                              r0=float(data[3]), desc=desc))
                    line = f.next()
            if line.strip() == 'ANGLES':
                line = f.next()
                while line.strip() != 'DIHEDRALS':
                    if line.strip() and line[0] != '!':
                        data, desc = line.split('!')[0], '!'.join(line.split('!')[1:]).strip()
                        data = data.split()
                        name = '{},{},{}'.format(data[0], data[1], data[2])
                        ff.angle_types.add(system.AngleType(name=name, tag=name, k=float(data[3]),
                                                                theta0=float(data[4]), desc=desc))
                    line = f.next()
            if line.strip() == 'DIHEDRALS':
                line = f.next()
                while line.strip() != 'IMPROPERS':
                    if line.strip() and line[0] != '!':
                        data, desc = line.split('!')[0], '!'.join(line.split('!')[1:]).strip()
                        data = data.split()
                        name = '{},{},{},{}'.format(data[0], data[1], data[2], data[3])
                        ff.dihedral_types.add(system.DihedralType(name=name, tag=name, k=float(data[4]),
                                                                      n=int(data[5]), d=cos(float(data[6])),
                                                                      desc=desc))
                    line = f.next()
            if line.strip() == 'IMPROPERS':
                line = f.next()
                while line.split() and line.split()[0] != 'NONBONDED':
                    if line.strip() and line[0] != '!':
                        data, desc = line.split('!')[0], '!'.join(line.split('!')[1:]).strip()
                        data = data.split()
                        name = '{},{},{},{}'.format(data[0], data[1], data[2], data[3])
                        ff.improper_types.add(system.ImproperType(name=name, tag=name, k=float(data[4]),
                                                                      x0=float(data[6]), desc=desc))
                    line = f.next()
            if line.split() and line.split()[0] == 'NONBONDED':
                line = f.next()
                line = f.next()
                line = f.next()
                while line.split() and line.split()[0] != 'HBOND':
                    if line.strip() and line[0] != '!':
                        data, desc = line.split('!')[0], '!'.join(line.split('!')[1:]).strip()
                        data = data.split()
                        name = '{}'.format(data[0])
                        if ff.particle_types.get(name):
                            pt = ff.particle_types.get(name)[0]
                            pt.epsilon = -1*float(data[2])
                            pt.sigma = float(data[3])*2/pow(2,1/6.)
                        else:
                            print name
                    line = f.next()
                    
    ff.write(out)

if __name__=='__main__':
    convert(sys.argv[1], sys.argv[2])
