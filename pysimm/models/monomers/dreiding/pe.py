from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk

def monomer():
    try:
        s = system.read_pubchem_smiles('CC')
    except:
        import os
        s = system.read_mol(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'CC.mol'))
    f = forcefield.Dreiding()
    s.apply_forcefield(f)
    
    c1 = s.particles[1]
    c2 = s.particles[2]
    c1.linker = 'head'
    c2.linker = 'tail'
    
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            s.particles.remove(pb.tag, update=False)
            break
        
    for b in c2.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c2 else b.b
            s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()
    
    s.pair_style = 'lj/cut'
    
    lmps.quick_min(s, min_style='fire')
    
    s.add_particle_bonding()
    
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