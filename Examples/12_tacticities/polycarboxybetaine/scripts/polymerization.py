from pysimm import lmps, system, forcefield
from pysimm.apps.random_walk import random_walk, random_walk_tacticity, check_tacticity

import numpy


# ----------> Here starts the main body of polymerisation code <------------
def cap_with_methyls(input_sst, ff):
    '''
        An utility method that implements capping of free ends of polymer chains with methyl 
        groups in all-atom forcefield representation
    '''
    # Let's cap the oligomer with the methyl (-CH3) group
    captypes = []
    for cpn in ['CG331', 'HGA3']:
        tmp = input_sst.particle_types.get(cpn)
        if tmp:
            cpt = tmp[0]
        else:
            cpt = ff.particle_types.get(cpn)[0].copy()
            input_sst.particle_types.add(cpt)
        captypes.append(cpt)

    for p in input_sst.particles:
        if p.linker is not None:
            if len(p.bonded_to) < 4:

                # assuming that the linker atom is sp3 hybridised C let's define the last non-occupied direction
                # of the tetrahedron
                dirctn = numpy.zeros(3)
                for p_ in p.bonded_to:
                    dirctn += numpy.array([p.x, p.y, p.z]) - numpy.array([p_.x, p_.y, p_.z])

                dirctn = dirctn / numpy.linalg.norm(dirctn)
                cap_c = system.Particle(x=p.x + 1.53 * dirctn[0], y=p.y + 1.53 * dirctn[1], z=p.z + 1.53 * dirctn[2],
                                        type=captypes[0])
                input_sst.add_particle_bonded_to(cap_c, p, f=ff)

                dir_h = numpy.array([1.0, 1.0, 1.0])
                dir_h[0] = -(dir_h[1] * dirctn[1] + dir_h[2] * dirctn[2]) / dirctn[0]
                dir_h = dir_h / numpy.linalg.norm(dir_h)

                dir_h2 = numpy.array([1.0, 1.0, -1.0])
                dir_h2[1] = (dirctn[2] / dirctn[0] - dir_h[2] / dir_h[0]) / (dirctn[1] / dirctn[0] - dir_h[1] / dir_h[0])
                dir_h2[0] = dirctn[2] / dirctn[0] - dirctn[1] * dir_h2[1] / dirctn[0]
                dir_h2 = dir_h2 / numpy.linalg.norm(dir_h2)

                stretch = 0.78
                input_sst.add_particle_bonded_to(system.Particle(x=cap_c.x + stretch * dirctn[0] + stretch * dir_h[0],
                                                                 y=cap_c.y + stretch * dirctn[1] + stretch * dir_h[1],
                                                                 z=cap_c.z + stretch * dirctn[2] + stretch * dir_h[2],
                                                                 type=captypes[1]), cap_c, f=ff)
                input_sst.add_particle_bonded_to(system.Particle(x=cap_c.x + stretch * dirctn[0] + stretch * dir_h2[0],
                                                                 y=cap_c.y + stretch * dirctn[1] + stretch * dir_h2[1],
                                                                 z=cap_c.z + stretch * dirctn[2] + stretch * dir_h2[2],
                                                                 type=captypes[1]), cap_c, f=ff)
                input_sst.add_particle_bonded_to(system.Particle(x=cap_c.x + stretch * dirctn[0] - stretch * dir_h2[0],
                                                                 y=cap_c.y + stretch * dirctn[1] - stretch * dir_h2[1],
                                                                 z=cap_c.z + stretch * dirctn[2] - stretch * dir_h2[2],
                                                                 type=captypes[1]), cap_c, f=ff)
    input_sst.objectify()
    input_sst.center(what='particles', at=[0.0, 0.0, 0.0], move_both=False)

    sim = lmps.Simulation(input_sst, log='capping_opt.log')
    sim.add_min(min_style='cg', name='min_cg', etol=1.0e-6, ftol=1.0e-6, maxiter=int(1e+6), maxeval=int(1e+7))
    sim.run()


# ----------> Here is the declaratin section -- majority of parameters <------------
# ---------->  to adjust for your monomer to polymerize it correctly   <------------

# defines whether free ends of oligomer chains will be capped with methyls or not
is_cap = True

# length of the polymer chain built in random_walk() polymerisation
chain_len = 10

# model of repititive unit for the polymer to built -- code assumes that the monomer is uncapped -- two atoms in
# the molecule which are called linkers are undercoordinated
data_path = '../../../../pysimm/models/monomers/topologies/'
monomer = system.read_pdb(data_path + 'cbma.pdb', str_file=data_path + 'cbma.str')

# mapping defines 'geometrically important' atoms in a monomer you use: 'head' and 'tail' are connection points for
# polymer chain growth; 'mirror' together with head and tail define plane which will mirrors the monomer for
# syndiotactic insertion. Best way to define the mirror atom is to assign it to the capping atom of the head, though
# the capping atom does not exist yet, it will be created before the random_walk_tacticity() run and will have next
# availible index (e.g. if there were 35 atoms in uncapped system it will have index #36)
lnkr_atoms = {'head': 1, 'tail': 2}

# list defines indexes of atoms for the tacticity analysis and marks 1st backbone atom, 2nd backbone atom, atom in
# the first side chain and atom in the second side chain (like methyl group)
tacticity_order = [1, 2, 3, 5]


# ----------> Here starts the main body of polymerisation code <------------
ff = forcefield.Charmm()

# we assign linkers by absolute indexes, so for your monomer they likely will be different
#  thus please assign them correct values
for nm in lnkr_atoms.keys():
    monomer.particles[lnkr_atoms[nm]].linker = nm

# Type forcefield particles types in the system automatically using Charmm-FF typer of PySIMM
monomer.apply_forcefield(ff, charges=None)

# In this example partial charges of the system are already set in the .str file and read into the system
# let's check whether they add up to 0, and add counterion if not, because some of the monomers in the library are charged.
monomer.set_charge()
print('\tRead monomer has charge of {:}q'.format(round(monomer.charge, 8)))
if abs(monomer.charge + 1) < 0.1:
    print('\tAdding SODIUM counterion to equilibrate the system chargewise')
    cntrion_tp = ff.particle_types.get('SOD')[0].copy()
    chrg = 1
elif abs(monomer.charge - 1) < 0.1:
    print('\tAdding CLORINE counterion to equilibrate the system chargewise')
    cntrion_tp = ff.particle_types.get('CLA')[0].copy()
    chrg = -1
else:
    cntrion_tp = None

if cntrion_tp:
    monomer.particle_types.add(cntrion_tp)
    monomer.particles.add(system.Particle(x=monomer.cog[0], y=monomer.cog[1], z=monomer.cog[2] + 5.0, type=cntrion_tp,
                                      charge=chrg, molecule=monomer.molecules[1]))

# -------------> Polymer construction and tacticity check <--------------
sngl_chain = random_walk(monomer, chain_len, forcefield=ff, density=0.01, print_to_screen='true', traj=False, unwrap=True)

if is_cap:
    cap_with_methyls(sngl_chain, ff)

# After polymerisation and possibly capping is done let's cleanup: remove all counterions and reshape the simulation
# box so that the chain is in the center of a cube with padding of 3 nm

sngl_chain.center(what='particles', at=[0.0, 0.0, 0.0], move_both=False)

if cntrion_tp:
    for p in sngl_chain.particles:
        if p.type.name == cntrion_tp.name:
            sngl_chain.particles.remove(p.tag, update=False)
            sngl_chain.molecules[p.molecule.tag].particles.remove(p.tag, update=False)
    sngl_chain.objectify()

bxSize = 30.0
for param in ['dx', 'dy', 'dz']:
    setattr(sngl_chain, param, bxSize)

sngl_chain.write_pdb('1.polymer.random_walk' + '.pdb')

# check tacticity of the simple polymer chain 
tacticity_stat = check_tacticity(sngl_chain, tacticity_order, len(monomer.particles))
print('\t Simple random_walk chain contains {:}/{:} meso- and {:}/{:} racemo- diads'.format(
                            tacticity_stat[1].count(True), chain_len, tacticity_stat[1].count(False), chain_len))

# >>>>>> Controlled tactisity <<<<<<<<<<<<
new_monomer = monomer.copy()

# random_walk_tacticity requires capped molecule, however, cap atoms will be deleted 
# those are baiscally dummy atoms, an can be of any type. Let's define them as a carbon
# backbone atoms, which will allow us not to change types of other carbon atoms in the backbone
captype = (new_monomer.particle_types.get('CG331') +
           new_monomer.particle_types.get('CG321') +
           new_monomer.particle_types.get('CG311'))[0]
# loop through the particles to add caps to linkers
for p in new_monomer.particles:
    if p.linker:
        # define and normalize directional vector for the capping atom
        tmp = numpy.array([sum([p.x - p_.x for p_ in p.bonded_to]), 
                           sum([p.y - p_.y for p_ in p.bonded_to]), 
                           sum([p.z - p_.z for p_ in p.bonded_to])])
        tmp = 1.54 * tmp / numpy.linalg.norm(tmp)
        # add new capping particle along the defined direction
        new_p = new_monomer.add_particle_bonded_to(system.Particle(x=p.x + tmp[0], y=p.y + tmp[1],
                                                                   z=p.z + tmp[2], type=captype), p, f=ff)
        # decorate particle with '****_cap' tag
        new_p.rnd_wlk_tag = p.linker + '_cap'
        if p.linker == 'head':
            setattr(new_p, 'linker', 'mirror')

new_monomer.objectified = False
new_monomer.objectify()

new_monomer.write_pdb('2.monomer_to_rndwlk-tacticity.pdb')

# first random walk **without** simulations
polymer = random_walk_tacticity(new_monomer, chain_len + 1, forcefield=ff, tacticity='syndiotactic', sim='no')
polymer.write_pdb('3.polymer_rndwlk-tacticity.no-sim.pdb')

tacticity_stat = check_tacticity(polymer, tacticity_order, len(monomer.particles))
print('\t Polymer chain built without MD contains {:}/{:} meso- and {:}/{:} racemo- diads'.format(
                            tacticity_stat[1].count(True), chain_len, tacticity_stat[1].count(False), chain_len))


# second, random walk **with** simulations
polymer = random_walk_tacticity(new_monomer, chain_len + 1, forcefield=ff, tacticity='syndiotactic')
polymer.write_pdb('4.polymer_rndwlk-tacticity.pdb')

tacticity_stat = check_tacticity(polymer, tacticity_order, len(monomer.particles))
print('\t Polymer chain built using MD contains {:}/{:} meso- and {:}/{:} racemo- diads'.format(
                            tacticity_stat[1].count(True), chain_len, tacticity_stat[1].count(False), chain_len))

