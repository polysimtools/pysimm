from pysimm import system
from pysimm import lmps
from pysimm.apps.random_walk import random_walk


def position_rule_pim1(prev, new):
    pass



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
                         extra_bonds=[1, 0], traj=False, unwrap=True)

mimi = lmps.Simulation(sngl_chain, print_to_screen=False, log='mimimization.log')
mimi.add(lmps.Init(cutoff=14.0, special_bonds='lj 0.0 0.0 0.0 coul 0.0 0.0 0.5', pair_modify={'mix': 'arithmetic'}))
mimi.add_min(min_style='sd', eol=1e-5, ftol=1e-5, maxiter=int(1e+3), maxeval=int(1e+5))
mimi.add_min(min_style='cg', eol=1e-5, ftol=1e-5, maxiter=int(1e+3), maxeval=int(1e+5))
mimi.run(np=4)

sngl_chain.write_lammps('test_ladder.lmps')

