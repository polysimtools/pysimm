from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk

def monomer():
    try:
        s = system.read_pubchem_smiles('CC(C)C(=O)OC')
    except:
        import os
        s = system.read_mol(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'CC(C)C(=O)OC.mol'))
    m = s.molecules[1]
    f = forcefield.Gaff2()
    
    s.apply_forcefield(f)
    
    c3 = s.particles[3]
    c4 = s.particles[4]
    
    for b in c3.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c3 else b.b
            s.particles.remove(pb.tag, update=False)
            break
        
    for b in c4.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c4 else b.b
            s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()

    s.set_box(padding=10)
    
    c3.linker = 'head'
    c4.linker = 'tail'
    
    lmps.quick_min(s, min_style='fire')
    
    s.add_particle_bonding()
    
    return s
    
def polymer_chain(length):
    mon = monomer()
    polym = random_walk(mon, length, forcefield=forcefield.Gaff2())
    return polym