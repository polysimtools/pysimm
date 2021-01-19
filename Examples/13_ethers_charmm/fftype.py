from pysimm import lmps
from pysimm import system
from pysimm import forcefield

sst = system.read_mol('ethylpropylether.mol')

for b in sst.bonds:
    b.order = 1

sst.apply_forcefield(forcefield.Charmm(), charges='gasteiger')
sst.set_charge()

sst.dim.dx = 25.0
sst.dim.dy = 25.0
sst.dim.dz = 25.0
sst.center('box', [0.0, 0.0, 0.0], True)

sst.write_lammps('to_sim.lmps')

sim = lmps.Simulation(sst, log='simulation.log', cutoff={'inner_lj': 10.0, 'lj': 12.0})

for nd_lj in sst.nondiag_lj_types:
    sim.add_custom('pair_coeff {} {} {}'.format(' '.join(map(str, nd_lj.atm_types)), nd_lj.epsilon, nd_lj.sigma))

sim.add_min(min_style = 'cg', name = 'grad_min',  etol = 1.0e-6, ftol = 1.0e-6)
sim.add_md(ensemble='nvt', timestep=0.1, length=int(2e+5))
sim.write_input()
sim.run()

sst.write_lammps('ethylpropylether.optimised.lmps')
