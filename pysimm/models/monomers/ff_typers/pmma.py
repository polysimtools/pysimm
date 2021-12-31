from pysimm import system, lmps
from pysimm.apps.random_walk import random_walk


def monomer(ff, is_capped=False):
    try:
        s = system.read_pubchem_smiles('CC(C)C(=O)OC')
    except:
        import os
        s = system.read_mol(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'topologies', 'CC(C)C(=O)OC.mol'))
    
    s.apply_forcefield(ff)
    
    c3 = s.particles[3]
    c4 = s.particles[4]

    if not is_capped:
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


def polymer_chain(length, ff):
    return random_walk(monomer(ff), length, forcefield=ff)