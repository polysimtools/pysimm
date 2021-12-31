from pysimm import system
from pysimm import lmps
from pysimm.apps.random_walk import random_walk
import numpy as np
import os


def my_rotation(direction, angle):

    Rth = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    Tz = np.eye(3)
    dist = np.linalg.norm(direction)
    Tz[0][0] = direction[2] / dist
    Tz[2][2] = direction[2] / dist
    Tz[0][2] = -np.sqrt(direction[0] ** 2 + direction[1] ** 2) / dist
    Tz[2][0] = np.sqrt(direction[0] ** 2 + direction[1] ** 2) / dist

    Txz = np.eye(3)
    pln_dist = np.sqrt(direction[0] ** 2 + direction[1] ** 2)
    Txz[0][0] = direction[0] / pln_dist
    Txz[1][1] = direction[0] / pln_dist
    Txz[0][1] = direction[1] / pln_dist
    Txz[1][0] = -direction[1] / pln_dist

    return np.linalg.inv(Txz) @ np.linalg.inv(Tz) @ Rth @ Tz @ Txz


def position_rule_pim1(new, prev):
    '''
    Redefines the position of the next PIM-1 monomer unit relative to subset of particles coordinates corresponding to
    the previous monomer unit. Initially both units have the same coordinates
    Args:
        new:
        prev:

    Returns:
        ordr: permutation of linkers identifiers that define hed-tail linker connection pairs
    '''
    new.add_particle_bonding()
    head_pos = np.zeros(3)
    tail_pos = np.zeros(3)
    bnd_lngth = 2

    for p in prev:
        if p.linker == 'head':
            head_pos += np.array([p.x, p.y, p.z])
        elif p.linker == 'tail':
            tail_pos += np.array([p.x, p.y, p.z])

    head_tags = [p.tag for p in new.particles if p.linker == 'head']
    tail_tags = [p.tag for p in new.particles if p.linker == 'tail']

    # 1. Displace the monomer by a vector <tails> - <heads> (<...> means averaging)
    displ = head_pos / len(head_pos) - tail_pos / len(tail_pos)
    displ_dir = displ / np.linalg.norm(displ)
    for p, p_ in zip(prev, new.particles):
        p_.x = p.x + displ[0] + bnd_lngth * displ_dir[0]
        p_.y = p.y + displ[1] + bnd_lngth * displ_dir[1]
        p_.z = p.z + displ[2] + bnd_lngth * displ_dir[2]

    # 2. With 1/2 probability perform rotation of the monomer to connect linkers in reverse order
    ordr = [1, 0]
    if np.random.random() > 0.5:
        ordr = [0, 1]
        rot_pnt_tag = 4 # id of the center-of-rotation atom
        rot_vect = [0, 0, 0] # normal of the rotation plane
        for tt in tail_tags:
            p = new.particles[tt]
            p_ = [pt for pt in new.particles[tt].bonded_to][0]
            rot_vect += np.array([p.x, p.y, p.z]) - np.array([p_.x, p_.y, p_.z])

        new.rotate(around=new.particles[rot_pnt_tag], rot_matrix=my_rotation(rot_vect, np.pi))

    # 3. rotate new monomer so that aromatic plane with heads is approximately parallel to the aromatic plane with
    # tails of previous monomer
    p = new.particles[tail_tags[0]]
    p_ = new.particles[tail_tags[1]]
    new.rotate(around=p, rot_matrix=my_rotation(np.array([p.x, p.y, p.z]) - np.array([p_.x, p_.y, p_.z]), -np.pi / 2))

    return ordr

lmps.FF_SETTINGS['trappe/amber']['dihedral_style'] = 'harmonic'
lmps.FF_SETTINGS['trappe/amber']['special_bonds'] = 'amber'

# Importing initial monomer system
mono = system.read_lammps('pim_unit.lmps', pair_style='lj/cut', bond_style='harmonic',
                          angle_style='harmonic', dihedral_style='harmonic')
mono.zero_charge()
mono.forcefield = 'trappe/amber'

# Decorating monomer with linker-type atoms
ptitles = [pt.name for pt in mono.particle_types]
tails = [24, 26]
for t in tails:
    mono.particles[t].linker = 'tail'
    mono.particles[t].type = mono.particle_types[ptitles.index(mono.particles[t].type.name[1:]) + 1]

heads = [28, 29]
for t in heads:
    mono.particles[t].linker = 'head'
    mono.particles[t].type = mono.particle_types[ptitles.index(mono.particles[t].type.name[1:]) + 1]

# Preliminary settings for building low density PIM-1 box with dispersing separate molecules over the whole box space
nchains = 10
ld_box_size = 65
chains = []
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
                                       'maxiter': int(1e+7), 'maxeval': int(1e+7)})
    chains.append(sngl_chain)
    # remove litter
    os.remove('polymer.xyz')
    os.remove('polymer.lmps')


pim1_ld = system.System()
for prp in ['ff_class', 'forcefield', 'pair_style', 'bond_style', 'angle_style', 'dihedral_style', 'improper_style']:
    setattr(pim1_ld, prp, getattr(mono, prp))

pim1_ld.dim.dx = ld_box_size
pim1_ld.dim.dy = ld_box_size
pim1_ld.dim.dz = ld_box_size

for idx in range(10):
    tmp = chains[idx].copy()
    tmp.shift_particles(*displ[2 * idx])
    pim1_ld.add(tmp, change_dim=False)

pim1_ld.wrap()
pim1_ld.set_density()
print('Density of the system is: {:} g/cc'.format(pim1_ld.density))
# pim1_ld = system.replicate(new_chains, [1] * len(new_chains), density=0.05, rand=False)

mimi = lmps.Simulation(pim1_ld, print_to_screen=False, log='mimimization.log')
mimi.add(lmps.Init(cutoff=14.0, special_bonds='amber', pair_modify={'mix': 'arithmetic'}))
mimi.add_min(min_style='sd', eol=1e-5, ftol=1e-5, maxiter=int(1e+5), maxeval=int(1e+6))
mimi.add_min(min_style='cg', eol=1e-5, ftol=1e-5, maxiter=int(1e+5), maxeval=int(1e+6))
mimi.run(np=6)

pim1_ld.write_lammps('pim1.smpl1.rndwlk.lmps')

