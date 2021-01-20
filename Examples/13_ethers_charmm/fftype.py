import numpy

from pysimm import lmps
from pysimm import system
from pysimm import forcefield

sst = system.read_mol('ethylpropylether.mol')

bxSize = 25.0

for b in sst.bonds:
    b.order = 1

ff = forcefield.Charmm()
sst.apply_forcefield(ff, charges='gasteiger')
sst.set_charge()

sst.dim.dx = bxSize
sst.dim.dy = bxSize
sst.dim.dz = bxSize
sst.center('box', [0.0, 0.0, 0.0], True)

solute = system.read_lammps('tipS3P.lmps', angle_style='charmm')
solute.set_mass()
# knowing volume of the box let's find out how many water molecules are needed:
ngrid = numpy.floor(((bxSize ** 3) * 0.6022 / solute.mass) ** (1.0 / 3.0))

rng = numpy.linspace(sst.dim.xlo, sst.dim.xhi, int(ngrid) + 1)
count = 0
for p in rng[:-1]:
    for q in rng[:-1]:
        for t in rng[:-1]:
            flags = []
            for prt in sst.particles:
                dist = numpy.linalg.norm(numpy.array([prt.x, prt.y, prt.z]) - numpy.array([p, q, t]))
                flags.append(dist > 1.7)

            if all(flags):
                count += 1
                tmp = solute.copy(dx=p, dy=q, dz=t)
                sst.add(tmp, change_dim=False, update_properties=False)

ff.assign_extra_ljtypes(sst)

# sst.write_lammps('to_sim.lmps')

# Create simulation and add directly add all nondiagonal LJ parameters from CHARMM (LAMMPS style)
sim = lmps.Simulation(sst, log='simulation.log', cutoff={'inner_lj': 10.0, 'lj': 12.0})
for nd_lj in sst.nondiag_lj_types:
    sim.add_custom('pair_coeff {} {} {}'.format(' '.join(map(str, nd_lj.atm_types)), nd_lj.epsilon, nd_lj.sigma))

# Define velocities
sim.add(lmps.Velocity(temperature=300.0, style='create'))

run_opts = {'name': 'main', 'pressure': {'iso': 'iso', 'damp': 500}, 'ensemble': 'npt',
            'timestep': 0.5, 'temperature': 300, 'run': int(6e+4)}
sim.add(lmps.OutputSettings(thermo={'args': ['step', 'time', 'temp', 'density', 'etotal', 'epair']}))

sim.add_custom('fix shck_fix all shake 0.001 20 0 b {:} a {:}\n'.format(sst.bond_types.get('H,O')[0].tag,
                                                                        sst.angle_types.get('H,O,H')[0].tag))
sim.add(lmps.MolecularDynamics(**run_opts))
sim.run()

sst.write_lammps('ethylpropylether.wtr_solution.lmps')
