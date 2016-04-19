#!/usr/bin/env python

##############################
#
#
#  forcefield.pcff Module
#
#
##############################

import sys
from xml.dom import minidom
from xml.etree import ElementTree as ET
from pysimm.system import *
from pysimm.forcefield import Pcff


def convert(file_, out):
    ff = Pcff(False)
    with file(file_) as f:
        for line in f:
            line = line.split()
            if len(line) > 0 and line[0] == '#atom_types':
                print('reading atom types...')
                for _ in range(7):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    desc = ' '.join(line[6:]).strip()
                    pt = ff.particle_types.get(line[2])
                    if pt:
                        pt = pt[0]
                        pt.mass = float(line[3])
                        pt.elem = line[4]
                        pt.desc = desc
                        pt.ref.append(int(line[1]))
                    else:
                        ff.particle_types.add(ParticleType(mass=float(line[3]),
                                                           elem=line[4],
                                                           name=line[2],
                                                           tag=line[2],
                                                           desc=desc,
                                                           ref=[line[1]]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#equivalence':
                print('reading equivalences')
                for _ in range(6):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    pt = ff.particle_types.get(line[2])
                    if pt:
                        pt = pt[0]
                        pt.eq_vwd = line[3]
                        pt.eq_bond = line[4]
                        pt.eq_angle = line[5]
                        pt.eq_dihedral = line[6]
                        pt.eq_improper = line[7]
                        pt.ref.append(int(line[1]))
                    else:
                        ff.particle_types.add(ParticleType(tag=line[2],
                                                           name=line[2],
                                                           eq_vwd=line[3],
                                                           eq_bond=line[4],
                                                           eq_angle=line[5],
                                                           eq_dihedral=line[6],
                                                           eq_improper=line[7],
                                                           ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#bond_increments':
                print('reading bond charge increments...')
                for _ in range(4):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    p1, p2, q1, q2 = line[2:6]
                    name = '%s,%s' % (p1, p2)
                    bt = ff.bond_types.get(name, order=True)
                    if bt:
                        bt = bt[0]
                        bt.q1 = q1
                        bt.q2 = q2
                        bt.ref.append(int(line[1]))
                    else:
                        ff.bond_types.add(BondType(tag=name,
                                                   name=name,
                                                   q1=q1,
                                                   q2=q2,
                                                   ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#quartic_bond':
                print('reading bonds...')
                for _ in range(6):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    name = ','.join(line[2:4])
                    r0, k2, k3, k4 = map(float, line[4:8])
                    ref = int(line[1])
                    bt = ff.bond_types.get(name)
                    if bt:
                        bt = bt[0]
                        bt.r0 = r0
                        bt.k2 = k2
                        bt.k3 = k3
                        bt.k4 = k4
                        bt.ref.append(ref)
                    else:
                        ff.bond_types.add(BondType(tag=name,
                                                   name=name,
                                                   r0=r0,
                                                   k2=k2,
                                                   k3=k3,
                                                   k4=k4,
                                                   ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#quartic_angle':
                print('reading angles...')
                for _ in range(7):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    name = ','.join(line[2:5])
                    theta0, k2, k3, k4 = map(float, line[5:9])
                    ref = int(line[1])
                    at = ff.angle_types.get(name)
                    if at:
                        at = at[0]
                        at.theta0 = theta0
                        at.k2 = k2
                        at.k3 = k3
                        at.k4 = k4
                        at.ref.append(ref)
                    else:
                        ff.angle_types.add(AngleType(tag=name,
                                                     name=name,
                                                     theta0=theta0,
                                                     k2=k2,
                                                     k3=k3,
                                                     k4=k4,
                                                     ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#torsion_3':
                print('reading dihedrals...')
                for _ in range(6):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    name = ','.join(line[2:6])
                    k1, phi1, k2, phi2, k3, phi3 = map(float, line[6:12])
                    ref = int(line[1])
                    dt = ff.dihedral_types.get(name)
                    if dt:
                        dt = dt[0]
                        dt.k1 = k1
                        dt.k2 = k2
                        dt.k3 = k3
                        dt.phi1 = phi1
                        dt.phi2 = phi2
                        dt.phi3 = phi3
                        dt.ref.append(ref)
                    else:
                        ff.dihedral_types.add(DihedralType(tag=name,
                                                           name=name,
                                                           k1=k1,
                                                           k2=k2,
                                                           k3=k3,
                                                           phi1=phi1,
                                                           phi2=phi2,
                                                           phi3=phi3,
                                                           ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#wilson_out_of_plane':
                print('reading impropers...')
                for _ in range(6):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    p2, p1, p3, p4 = line[2:6]
                    name = ','.join([p1, p2, p3, p4])
                    k, x0 = map(float, line[6:8])
                    ref = int(line[1])
                    it = ff.improper_types.get(name)
                    if it:
                        it = it[0]
                        it.k = k
                        it.x0 = x0
                        it.ref.append(ref)
                    else:
                        ff.improper_types.add(ImproperType(tag=name,
                                                           name=name,
                                                           k=k,
                                                           x0=x0,
                                                           ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#nonbond(9-6)':
                print('reading nondonded terms...')
                for _ in range(13):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    name = line[2]
                    sigma = float(line[3])
                    epsilon = float(line[4])
                    pt = ff.particle_types.get(name)
                    if pt:
                        pt = pt[0]
                        pt.sigma = sigma
                        pt.epsilon = epsilon
                        pt.ref.append(int(line[1]))
                    else:
                        ff.particle_types.add(ParticleType(tag=name,
                                                           name=name,
                                                           sigma=sigma,
                                                           epsilon=epsilon,
                                                           ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#bond-bond':
                print('reading bond-bond angles...')
                for _ in range(6):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    p1, p2, p3 = line[2:5]
                    name = ','.join([p1, p2, p3])
                    m = float(line[5])
                    at = ff.angle_types.get(name)
                    if at:
                        at = at[0]
                        at.m = m
                        at.ref.append(int(line[1]))
                    else:
                        ff.angle_types.add(AngleType(tag=name,
                                                     name=name,
                                                     m=m,
                                                     ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#bond-bond_1_3':
                print('reading bond-bond-1-3 dihedrals...')
                for _ in range(6):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    p1, p2, p3, p4 = line[2:6]
                    name = ','.join([p1, p2, p3, p4])
                    n = float(line[6])
                    dt = ff.dihedral_types.get(name)
                    if dt:
                        dt = dt[0]
                        dt.n = n
                        dt.ref.append(int(line[1]))
                    else:
                        ff.dihedral_types.add(DihedralType(tag=name,
                                                           name=name,
                                                           n=n,
                                                           ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#bond-angle':
                print('reading bond-angle angles...')
                for _ in range(6):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    p1, p2, p3 = line[2:5]
                    name = ','.join([p1, p2, p3])
                    n1 = float(line[5])
                    n2 = float(line[6]) if len(line) == 7 else float(line[5])
                    at = ff.angle_types.get(name)
                    if at:
                        at = at[0]
                        at.n1 = n1
                        at.n2 = n2
                        at.ref.append(int(line[1]))
                    else:
                        ff.angle_types.add(AngleType(tag=name,
                                                     name=name,
                                                     n1=n1,
                                                     n2=n2,
                                                     ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#angle-angle':
                print('reading angle-angle impropers...')
                for _ in range(6):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    p1, p2, p3, p4 = line[2:6]
                    n1 = ','.join([p1, p2, p3, p4])
                    n2 = ','.join([p3, p2, p1, p4])
                    n3 = ','.join([p1, p2, p4, p3])
                    m = float(line[6])
                    it1 = ff.improper_types.get(n1)
                    if it1:
                        it1[0].m1 = m
                        it1[0].ref.append(int(line[1]))
                    else:
                        ff.improper_types.add(ImproperType(name=n1,
                                                           tag=n1,
                                                           m1=m,
                                                           ref=[int(line[1])]))
                    it2 = ff.improper_types.get(n2)
                    if it2:
                        it2[0].m2 = m
                        it2[0].ref.append(int(line[1]))
                    else:
                        ff.improper_types.add(ImproperType(tag=n2,
                                                           name=n2,
                                                           m2=m,
                                                           ref=[int(line[1])]))
                    it3 = ff.improper_types.get(n3)
                    if it3:
                        it3[0].m3 = m
                        it3[0].ref.append(int(line[1]))
                    else:
                        ff.improper_types.add(ImproperType(tag=n3,
                                                           name=n3,
                                                           m3=m,
                                                           ref=[int(line[1])]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#end_bond-torsion_3':
                print('reading end-bond-torsion dihedrals...')
                for _ in range(8):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    p1, p2, p3, p4 = line[2:6]
                    name = ','.join([p1, p2, p3, p4])
                    ref = int(line[1])
                    b1, b2, b3 = map(float, line[6:9])
                    if len(line) == 12:
                        c1, c2, c3 = map(float, line[9:12])
                    else:
                        c1, c2, c3 = b1, b2, b3
                    dt = ff.dihedral_types.get(name)
                    if dt:
                        dt = dt[0]
                        d_ = {'b1': b1, 'b2': b2, 'b3': b3,
                              'c1': c1, 'c2': c2, 'c3': c3}
                        dt.set(**d_)
                        dt.ref.append(ref)
                    else:
                        ff.dihedral_types.add(DihedralType(tag=name,
                                                           name=name,
                                                           b1=b1,
                                                           b2=b2,
                                                           b3=b3,
                                                           c1=c1,
                                                           c2=c2,
                                                           c3=c3,
                                                           ref=[ref]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#middle_bond-torsion_3':
                print('reading middle-bond-torsion dihedrals...')
                for _ in range(7):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    p1, p2, p3, p4 = line[2:6]
                    name = ','.join([p1, p2, p3, p4])
                    ref = int(line[1])
                    a1, a2, a3 = map(float, line[6:9])
                    dt = ff.dihedral_types.get(name)
                    if dt:
                        dt = dt[0]
                        dt.a1 = a1
                        dt.a2 = a2
                        dt.a3 = a3
                        dt.ref.append(ref)
                    else:
                        ff.dihedral_types.add(DihedralType(tag=name,
                                                           name=name,
                                                           a1=a1,
                                                           a2=a2,
                                                           a3=a3,
                                                           ref=[ref]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#angle-torsion_3':
                print('reading angle-torsion dihedrals...')
                for _ in range(9):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    p1, p2, p3, p4 = line[2:6]
                    name = ','.join([p1, p2, p3, p4])
                    ref = int(line[1])
                    d1, d2, d3 = map(float, line[6:9])
                    if len(line) == 12:
                        e1, e2, e3 = map(float, line[9:12])
                    else:
                        e1, e2, e3 = d1, d2, d3
                    dt = ff.dihedral_types.get(name)
                    if dt:
                        dt = dt[0]
                        d_ = {'d1': d1, 'd2': d2, 'd3': d3,
                              'e1': e1, 'e2': e2, 'e3': e3}
                        dt.set(**d_)
                        dt.ref.append(ref)
                    else:
                        ff.dihedral_types.add(DihedralType(tag=name,
                                                           name=name,
                                                           d1=d1,
                                                           d2=d2,
                                                           d3=d3,
                                                           e1=e1,
                                                           e2=e2,
                                                           e3=e3,
                                                           ref=[ref]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

            if len(line) > 0 and line[0] == '#angle-angle-torsion_1':
                print('reading angle-angle-torsion dihedrals...')
                for _ in range(6):
                    line = f.next().split()
                while line:
                    if line[0] == '#':
                        break
                    p1, p2, p3, p4 = line[2:6]
                    name = ','.join([p1, p2, p3, p4])
                    ref = int(line[1])
                    m = float(line[6])
                    dt = ff.dihedral_types.get(name)
                    if dt:
                        dt = dt[0]
                        dt.m = m
                        dt.ref.append(ref)
                    else:
                        ff.dihedral_types.add(DihedralType(tag=name,
                                                           name=name,
                                                           m=m,
                                                           ref=[ref]))
                    line = f.next()
                    if len(line.strip()) == 0:
                        line = f.next()
                    line = line.split()

    print('setting extra cross term parameters r1 and r2 for angle types')
    for at in ff.angle_types:
        p1, p2, p3 = at.name.split(',')
        p1_ = ff.particle_types.get(p1)
        if p1_:
            p1 = p1_[0].eq_bond or p1_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p1)
            continue
        p2_ = ff.particle_types.get(p2)
        if p2_:
            p2 = p2_[0].eq_bond or p2_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p2)
            continue
        p3_ = ff.particle_types.get(p3)
        if p3_:
            p3 = p3_[0].eq_bond or p3_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p3)
            continue

        b1 = ff.bond_types.get(','.join([p1, p2]))
        if b1:
            at.r1 = b1[0].r0
        else:
            print('cannot find r0 for BondType %s...skipping'
                  % ','.join([p1, p2]))
        b2 = ff.bond_types.get(','.join([p2, p3]))
        if b2:
            at.r2 = b2[0].r0
        else:
            print('cannot find r0 for BondType %s...skipping' %
                  ','.join([p2, p3]))

    print('setting extra cross term parameters r1, r2, r3, '
          'theta1, theta2 for dihedral types')
    for dt in ff.dihedral_types:
        p1, p2, p3, p4 = dt.name.split(',')
        p1_ = ff.particle_types.get(p1)
        if p1_:
            p1b = p1_[0].eq_bond or p1_[0].name
            p1a = p1_[0].eq_angle or p1_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p1)
            continue
        p2_ = ff.particle_types.get(p2)
        if p2_:
            p2b = p2_[0].eq_bond or p2_[0].name
            p2a = p2_[0].eq_angle or p2_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p2)
            continue
        p3_ = ff.particle_types.get(p3)
        if p3_:
            p3b = p3_[0].eq_bond or p3_[0].name
            p3a = p3_[0].eq_angle or p3_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p3)
            continue
        p4_ = ff.particle_types.get(p4)
        if p4_:
            p4b = p4_[0].eq_bond or p4_[0].name
            p4a = p4_[0].eq_angle or p4_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p4)
            continue

        b1 = ff.bond_types.get(','.join([p1b, p2b]))
        if b1:
            dt.r1 = b1[0].r0
        else:
            print('cannot find r0 for BondType %s...skipping'
                  % ','.join([p1b, p2b]))
        b2 = ff.bond_types.get(','.join([p2b, p3b]))
        if b2:
            dt.r2 = b2[0].r0
        else:
            print('cannot find r0 for BondType %s...skipping'
                  % ','.join([p2b, p3b]))
        b3 = ff.bond_types.get(','.join([p3b, p4b]))
        if b3:
            dt.r3 = b3[0].r0
        else:
            print('cannot find r0 for AngleType %s...skipping'
                  % ','.join([p3b, p4b]))
        a1 = ff.angle_types.get(','.join([p1a, p2a, p3a]))
        if a1:
            dt.theta1 = a1[0].theta0
        else:
            print('cannot find r0 for AngleType %s...skipping'
                  % ','.join([p1a, p2a, p3a]))
        a2 = ff.angle_types.get(','.join([p2a, p3a, p4a]))
        if a2:
            dt.theta2 = a2[0].theta0
        else:
            print('cannot find r0 for AngleType %s...skipping'
                  % ','.join([p2a, p3a, p4a]))

    print('setting extra cross term parameters theta1, theta2, theta3 for '
          'improper types')
    for it in ff.improper_types:
        p2, p1, p3, p4 = it.name.split(',')
        p1_ = ff.particle_types.get(p1)
        if p1_:
            p1_ = p1_[0].eq_improper or p1_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p1)
            continue
        p2_ = ff.particle_types.get(p2)
        if p2_:
            p2_ = p2_[0].eq_improper or p2_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p2)
            continue
        p3_ = ff.particle_types.get(p3)
        if p3_:
            p3_ = p3_[0].eq_improper or p3_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p3)
            continue
        p4_ = ff.particle_types.get(p4)
        if p4_:
            p4_ = p4_[0].eq_improper or p4_[0].name
        else:
            print('cannot find ParticleType %s...skipping' % p4)
            continue

        a1 = ff.angle_types.get(','.join([p1_, p2_, p3_]))
        if a1:
            it.theta1 = a1[0].theta0
        else:
            print('cannot find theta0 for AngleType %s...skipping'
                  % ','.join([p1_, p2_, p3_]))

        a2 = ff.angle_types.get(','.join([p1_, p2_, p4_]))
        if a2:
            it.theta2 = a2[0].theta0
        else:
            print('cannot find theta0 for AngleType %s...skipping'
                  % ','.join([p1_, p2_, p4_]))

        a3 = ff.angle_types.get(','.join([p3_, p2_, p4_]))
        if a3:
            it.theta3 = a3[0].theta0
        else:
            print('cannot find theta0 for AngleType %s...skipping'
                  % ','.join([p3_, p2_, p4_]))



    ff.write(out)


if __name__ == '__main__':
    convert(sys.argv[1], sys.argv[2])
