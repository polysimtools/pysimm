from pysimm import system, lmps, forcefield

def run(test=False):
    # create empty system
    print('Example progress: Creating an empty system...')
    s = system.System()
    
    # create new molecule in our system
    print('Example progress: Adding an empty molecule container to our system...')
    m = s.molecules.add(system.Molecule())
    
    # retrieve Dreiding parameters
    print('Example progress: Retrieving Dreiding force field parameters...')
    f = forcefield.Dreiding()
    s.forcefield = f.name
    
    # get a copy of the C_ particle type object from Dreiding
    # get method returns a list, we need the first element
    dreiding_C_ = s.particle_types.add(f.particle_types.get('C_3')[0].copy())
    
    # get H_ particle type object from Dreiding
    dreiding_H_ = s.particle_types.add(f.particle_types.get('H_')[0].copy())
    
    # we'll first make the carbon atom at the origin
    # we'll include gasteiger charges later
    print('Example progress: Adding carbon atom at origin...')
    c1 = s.particles.add(system.Particle(type=dreiding_C_, x=0, y=0, z=0, charge=0, molecule=m))
    
    # now we'll add 4 hydrogen atoms bonded to our carbon atom
    # these atoms will be placed randomly 1.5 angstroms from the carbon atom
    # we'll optimize the structure using LAMMPS afterwords
    # we supply the Dreiding forcefield object so that bond and angle types can be added as well
    print('Example progress: Adding 4 hydrogen atoms at random positions bonded to the carbon atom...')
    h1 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
    h2 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
    h3 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
    h4 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
    
    # let's add gasteiger charges
    print('Example progress: Deriving Gasteiger charges...')
    s.apply_charges(f, charges='gasteiger')
    
    # right now there is no simulation box defined
    # we'll define a box surrounding our methane molecule with a 10 angstrom padding
    print('Example progress: Constructing Simulation box surrounding our new molecule...')
    s.set_box(padding=10)
    
    # before we optimize our structure, LAMMPS needs to know what type of 
    # pair, bond, and angle interactions we are using
    # these are determined by the forcefield being used
    s.pair_style='lj/cut'
    s.bond_style='harmonic'
    s.angle_style='harmonic'
    
    # we'll perform energy minimization using the fire algorithm in LAMMPS
    print('Example progress: Optimizing structure using LAMMPS...')
    lmps.quick_min(s, min_style='fire', name='fire_min', etol=1e-10, ftol=1e-10)
    
    # write xyz, YAML, LAMMPS data, and chemdoodle json files
    print('Example progress: Saving structure to files...')
    s.write_xyz('methane.xyz')
    s.write_yaml('methane.yaml')
    s.write_lammps('methane.lmps')
    s.write_chemdoodle_json('methane.json')
    
    print('Example progress: Complete!')
    
if __name__ == '__main__':
    run()
