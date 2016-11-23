from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk

def monomer():
    s = system.System()
    m = s.molecules.add(system.Molecule())
    f = forcefield.Gaff2()
    
    s.pair_style = f.pair_style
    s.bond_style = f.bond_style
    s.angle_style = f.angle_style
    s.dihedral_style = f.dihedral_style
    s.improper_style = f.improper_style
    
    gaff_c3 = s.particle_types.add(f.particle_types.get('c3')[0].copy())
    c1 = s.particles.add(system.Particle(type=gaff_c3, x=0, y=0, z=0, charge=0, molecule=m))
    c2 = s.add_particle_bonded_to(system.Particle(type=gaff_c3, charge=0, molecule=m), c1, f)
    
    gaff_hc = s.particle_types.add(f.particle_types.get('hc')[0].copy())
    
    h1 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c1, f)
    h2 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c1, f)
    h3 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c2, f)
    h4 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c2, f)

    s.set_box(padding=10)
    
    c1.linker = 'head'
    c2.linker = 'tail'
    
    lmps.quick_min(s, min_style='fire')
    
    return s
    
def polymer_chain(length):
    mon = monomer()
    polym = random_walk(mon, length, forcefield=forcefield.Dreiding())
    return polym