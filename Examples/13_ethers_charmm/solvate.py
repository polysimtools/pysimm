import numpy
from pysimm import lmps
from pysimm import system
from pysimm import forcefield

sst = system.read_mol('ethylpropylether.mol')

# decorate the system read from a simple .mol file and type it with CHARMM FF parameters and simple gasteiger charges
for b in sst.bonds:
    b.order = 1
ff = forcefield.Charmm()
sst.apply_forcefield(ff, charges='gasteiger')
sst.set_charge()

# also expand the system to a minimal size at which the periodic image convention is held
bxSize = 25.0
sst.dim.dx = bxSize
sst.dim.dy = bxSize
sst.dim.dz = bxSize
sst.center('box', [0.0, 0.0, 0.0], True)

# read a water molecule from .lmps file that already contains all details
solvnt = system.read_lammps('tipS3P.lmps', angle_style='charmm')

# knowing volume of the box let's find out how many water molecules are needed to feet at each side of a cubic box:
solvnt.set_mass()
ngrid = numpy.floor(((bxSize ** 3) * 0.6022 / solvnt.mass) ** (1.0 / 3.0))

# put water molecules in the nodes of a regular grid
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
                tmp = solvnt.copy(dx=p, dy=q, dz=t)
                sst.add(tmp, change_dim=False, update_properties=False)

# because 2 systems have been combined let's reassign non-diagonal LJ interactions
ff.assign_extra_ljtypes(sst)

sst.write_lammps('to_sim.lmps')

# Create simulation and directly add all nondiagonal LJ parameters to the run file
# (that is how LAMMPS operates)
sim = lmps.Simulation(sst, log='simulation.log', cutoff={'inner_lj': 10.0, 'lj': 12.0})
if sst.nondiag_lj_types:
    for nd_lj in sst.nondiag_lj_types:
        sim.add_custom('pair_coeff {} {} {}'.format(' '.join(map(str, nd_lj.atm_types)), nd_lj.epsilon, nd_lj.sigma))

# define velocities and output settings
sim.add(lmps.Velocity(temperature=300.0, style='create'))
sim.add(lmps.OutputSettings(thermo={'args': ['step', 'time', 'temp', 'density', 'etotal', 'epair']}))

# setup shake for H-O bond and H-O-H angle for the whole simulation system
sim.add_custom('fix shck_fix all shake 0.001 40 0 b {:} a {:}\n'.format(sst.bond_types.get('H,O')[0].tag,
                                                                        sst.angle_types.get('H,O,H')[0].tag))
sim.add(lmps.MolecularDynamics(name='main',
                               pressure={'iso': 'iso', 'damp': 100},
                               ensemble='npt',
                               timestep=1,
                               temperature=300.0,
                               run=int(5e+4)))
sim.run()

sst.write_lammps('eth_propether.wtr_solution.lmps')

