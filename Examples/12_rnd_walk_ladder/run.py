from pysimm import system
from pysimm import lmps
from pysimm.apps.random_walk import random_walk
import numpy as np
import os

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
lmps.FF_SETTINGS['trappe/amber']['special_bonds'] = 'amber'

ptittles = [pt.name for pt in mono.particle_types]

tails = [24, 26]
for t in tails:
    mono.particles[t].linker = 'tail'
    mono.particles[t].type = mono.particle_types[ptittles.index(mono.particles[t].type.name[1:]) + 1]

heads = [28, 29]
for t in heads:
    mono.particles[t].linker = 'head'
    mono.particles[t].type = mono.particle_types[ptittles.index(mono.particles[t].type.name[1:]) + 1]

nchains = 10
chains = []
mono.zero_charge()

x = np.linspace(-32, 32, 3)
X, Y, Z = np.meshgrid(x, x, x)
displ = [[x, y, z] for x, y, z in zip(np.ndarray.flatten(X), np.ndarray.flatten(Y), np.ndarray.flatten(Z))]

for ch in range(1, nchains + 1):

    mloc = mono.copy()
    mloc.dim.translate(*displ[2 * ch])
    mloc.shift_particles(*displ[2 * ch])
    
    sngl_chain = random_walk(mloc, 10, density=0.01,
                             geometry_rule=position_rule_pim1, extra_bonds=[1, 0], traj=False, unwrap=True,
                             settings={'np': 6, 'length': 25000, 'min_style': 'cg', 'etol': 1e-5, 'ftol': 1e-5,
                                       'maxiter': int(1e+5), 'maxeval': int(1e+6)})

    chains.append(sngl_chain)
    # remove litter
    os.remove('polymer.xyz')
    os.remove('polymer.lmps')



pim1_ld = system.System()
pim1_ld.ff_class = sngl_chain.ff_class
pim1_ld.forcefield = sngl_chain.forcefield
pim1_ld.pair_style = sngl_chain.pair_style
pim1_ld.bond_style = sngl_chain.bond_style
pim1_ld.angle_style = sngl_chain.angle_style
pim1_ld.dihedral_style = sngl_chain.dihedral_style
pim1_ld.improper_style = sngl_chain.improper_style

pim1_ld.dim.dx = 65
pim1_ld.dim.dy = 65
pim1_ld.dim.dz = 65

for idx in range(10):
    tmp = chains[idx].copy()
    tmp.shift_particles(*displ[2 * idx])
    pim1_ld.add(tmp, change_dim=False)


pim1_ld.wrap()
pim1_ld.set_density()
print(pim1_ld.density)
# pim1_ld = system.replicate(new_chains, [1] * len(new_chains), density=0.05, rand=False)


mimi = lmps.Simulation(pim1_ld, print_to_screen=False, log='mimimization.log')
mimi.add(lmps.Init(cutoff=14.0, special_bonds='amber', pair_modify={'mix': 'arithmetic'}))
mimi.add_min(min_style='sd', eol=1e-5, ftol=1e-5, maxiter=int(1e+5), maxeval=int(1e+6))
mimi.add_min(min_style='cg', eol=1e-5, ftol=1e-5, maxiter=int(1e+5), maxeval=int(1e+6))
mimi.run(np=6)


pim1_ld.write_lammps('pim1.smpl5.rndwlk.lmps')

