from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk

def monomer():
    s = system.System()
    m = s.molecules.add(system.Molecule())
    f = forcefield.Dreiding()
    
    s.pair_style = f.pair_style
    s.bond_style = f.bond_style
    s.angle_style = f.angle_style
    s.dihedral_style = f.dihedral_style
    s.improper_style = f.improper_style
    
    dreiding_C_3 = s.particle_types.add(f.particle_types.get('C_3')[0].copy())
    c1 = s.particles.add(system.Particle(type=dreiding_C_3, x=0, y=0, z=0, charge=0, molecule=m))
    c2 = s.add_particle_bonded_to(system.Particle(type=dreiding_C_3, charge=0, molecule=m), c1, f)
    
    dreiding_H_ = s.particle_types.add(f.particle_types.get('H_')[0].copy())
    
    h1 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
    h2 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
    h3 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c2, f)
    h4 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c2, f)
    
    s.set_box(padding=10)
    
    c1.linker = 'head'
    c2.linker = 'tail'
    
    s.pair_style = 'lj'
    
    lmps.quick_min(s, min_style='fire')
    
    return s
    
def polymer_chain(length):
    mon = monomer()
    polym = random_walk(mon, length, forcefield=forcefield.Dreiding())
    return polym
    
def polymer_system(chains=10, mn=1000, pdi=1, density=0.3):
    if pdi != 1:
        print('disperse molecular weight distributions not supported yet')
        return
    
    mon = monomer()
    
    chain_length = int(mn/mon.mass)
    
    polym = random_walk(mon, chain_length, density=density/chains, forcefield=forcefield.Dreiding())
    
    for chain in range(chains-1):
        polym = random_walk(mon, chain_length, s_=polym, density=None, forcefield=forcefield.Dreiding())
    
    return polym