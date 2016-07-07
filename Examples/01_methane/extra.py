from pysimm import system, lmps, forcefield, gasteiger

# create empty system
s = system.System()

# create new molecule in our system
m = s.molecules.add(system.Molecule())

# retrieve GAFF parameters
f = forcefield.Gaff()

# get a copy of the c3 particle type object from GAFF
# get method returns a list, we need the first element
gaff_c3 = f.particle_types.get('c3')[0].copy()

# now we need to add this particle type to our system
s.particle_types.add(gaff_c3)

# we'll first make the carbon atom at the origin
# we'll include gasteiger charges later
c1 = s.particles.add(system.Particle(type=gaff_c3, x=0, y=0, z=0, charge=0, molecule=m))

# get hc particle type object from GAFF
gaff_hc = f.particle_types.get('hc')[0].copy()

# add hc particle type to system
s.particle_types.add(gaff_hc)

# now we'll add 4 hydrogen atoms bonded to our carbon atom
# these atoms will be placed randomly 1.5 angstroms from the carbon atom
# we'll optimize the structure using LAMMPS afterwords
# we supply the GAFF forcefield object so that bond and angle types can be added as well
h1 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c1, f)
h2 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c1, f)
h3 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c1, f)

c2 = s.add_particle_bonded_to(system.Particle(type=gaff_c3, charge=0, molecule=m), c1, f)

h4 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c2, f)
h5 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c2, f)
h6 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c2, f)

# right now there is no simulation box defined
# we'll define a box surrounding our methane molecule with a 10 angstrom padding
s.set_box(padding=10)

# before we optimize our structure, LAMMPS needs to know what type of 
# pair, bond, and angle interactions we are using
# these are determined by the forcefield being used
s.pair_style='lj'
s.bond_style='harmonic'
s.angle_style='harmonic'
s.dihedral_style='harmonic'

# we'll do a two step energy minimization using LAMMPS
# step 1: steepest decent algorithm
# step 2: conjugent gradient algorithm
lmps.quick_min(s, min_style='sd', print_to_screen=True)
lmps.quick_min(s, min_style='cg', print_to_screen=True)

# let's add gasteiger charges
gasteiger.set_charges(s)

# write a YAML file and a LAMMPS data file
s.write_yaml('ethane.yaml')
s.write_lammps('ethane.lmps')

# take a look at our new molecule in VMD
s.viz()
