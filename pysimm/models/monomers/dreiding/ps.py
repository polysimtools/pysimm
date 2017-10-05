from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk

def monomer():
    try:
        s = system.read_pubchem_smiles('CCc1=cc=cc=c1')
    except:
        import os
        s = system.read_mol(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'CCc1=cc=cc=c1.mol'))
    m = s.molecules[1]
    f = forcefield.Dreiding()
    
    for b in s.bonds:
        if b.a.bonds.count == 3 and b.b.bonds.count == 3:
            b.order = 4
    
    s.apply_forcefield(f)
    
    c1 = s.particles[1]
    c5 = s.particles[5]
    
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            s.particles.remove(pb.tag, update=False)
            break
        
    for b in c5.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c5 else b.b
            s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()

    s.set_box(padding=10)
    
    c1.linker = 'head'
    c5.linker = 'tail'
    
    lmps.quick_min(s, min_style='fire')
    
    s.add_particle_bonding()
    
    return s
    
def polymer_chain(length):
    mon = monomer()
    polym = random_walk(mon, length, forcefield=forcefield.Dreiding())
    return polym
    
    