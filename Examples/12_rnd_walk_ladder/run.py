from pysimm import system
from pysimm import lmps
from pysimm.apps.random_walk import random_walk
import numpy as np


def position_rule_pim1(new, prev):

    head_pos = np.zeros(3)
    tail_pos = np.zeros(3)
    hcount = 0
    tcount = 0
    bnd_lngth = 1.8

    for p in prev:
        if p.linker == 'head':
            hcount += 1
            head_pos += np.array([p.x, p.y, p.z])
        elif p.linker == 'tail':
            tcount += 1
            tail_pos += np.array([p.x, p.y, p.z])

    displ = head_pos / hcount - tail_pos / tcount
    displ_dir = displ / np.linalg.norm(displ)
    for p, p_ in zip(prev, new.particles):
        p_.x = p.x + displ[0] + bnd_lngth * displ_dir[0]
        p_.y = p.y + displ[1] + bnd_lngth * displ_dir[1]
        p_.z = p.z + displ[2] + bnd_lngth * displ_dir[2]

    is_rotation = np.random.random() > 0.5
    # is_rotation = True
    ordr = [1, 0]
    if is_rotation:
        ordr = [0, 1]
        # Creating an arbitrary rotation of added monomer around the common backbone direction
        rot_pnt_tag = 4

        theta = np.pi
        Rth = np.eye(3)
        Rth[0][0] = np.cos(theta)
        Rth[1][1] = np.cos(theta)
        Rth[0][1] = -np.sin(theta)
        Rth[1][0] = np.sin(theta)

        Tz = np.eye(3)

        new.add_particle_bonding()
        rot_vect = [0, 0, 0]
        for p in new.particles:
            if p.linker == 'tail':
                for p_ in p.bonded_to:
                    rot_vect += np.array([p.x, p.y, p.z]) - np.array([p_.x, p_.y, p_.z])

        dist = np.linalg.norm(rot_vect)
        Tz[0][0] = rot_vect[2] / dist
        Tz[2][2] = rot_vect[2] / dist
        Tz[0][2] = -np.sqrt(rot_vect[0] ** 2 + rot_vect[1] ** 2) / dist
        Tz[2][0] = np.sqrt(rot_vect[0] ** 2 + rot_vect[1] ** 2) / dist

        Txz = np.eye(3)
        pln_dist = np.sqrt(rot_vect[0] ** 2 + rot_vect[1] ** 2)
        Txz[0][0] = rot_vect[0] / pln_dist
        Txz[1][1] = rot_vect[0] / pln_dist
        Txz[0][1] = rot_vect[1] / pln_dist
        Txz[1][0] = -rot_vect[1] / pln_dist

        new.rotate(around=new.particles[rot_pnt_tag], rot_matrix=np.linalg.inv(Txz) @ np.linalg.inv(Tz) @ Rth @ Tz @ Txz)
    return ordr


mono = system.read_lammps('pim_unit.lmps', pair_style='lj/cut', bond_style='harmonic',
                          angle_style='harmonic', dihedral_style='harmonic')
mono.forcefield = 'trappe/amber'
lmps.FF_SETTINGS['trappe/amber']['dihedral_style'] = 'harmonic'


ptittles = [pt.name for pt in mono.particle_types]

tails = [24, 26]
for t in tails:
    mono.particles[t].linker = 'tail'
    mono.particles[t].type = mono.particle_types[ptittles.index(mono.particles[t].type.name[1:]) + 1]

heads = [28, 29]
for t in heads:
    mono.particles[t].linker = 'head'
    mono.particles[t].type = mono.particle_types[ptittles.index(mono.particles[t].type.name[1:]) + 1]

sngl_chain = random_walk(mono, 50, density=0.01,
                         geometry_rule=position_rule_pim1, extra_bonds=[1, 0], traj=False, unwrap=True,
                         settings={'length': 10000, 'min_style': 'cg', 'etol': 1e-5, 'ftol': 1e-5, 'maxiter': int(1e+4), 'maxeval': int(1e+5)})


mimi = lmps.Simulation(sngl_chain, print_to_screen=False, log='mimimization.log')
mimi.add(lmps.Init(cutoff=14.0, special_bonds='lj 0.0 0.0 0.0 coul 0.0 0.0 0.0', pair_modify={'mix': 'arithmetic'}))
mimi.add_min(min_style='sd', eol=1e-5, ftol=1e-5, maxiter=int(1e+4), maxeval=int(1e+5))
mimi.add_min(min_style='cg', eol=1e-5, ftol=1e-5, maxiter=int(1e+4), maxeval=int(1e+5))
mimi.run(np=4)

sngl_chain.write_lammps('test_ladder.lmps')

