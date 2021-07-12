from pysimm import system, lmps
from pysimm.apps.random_walk import random_walk


def monomer(ff, is_capped=False):

    try:
        s = system.read_pubchem_smiles('CCc1=cc=cc=c1')
    except:
        import os
        s = system.read_mol(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'topologies', 'CCc1=cc=cc=c1.mol'))

    for b in s.bonds:
        if b.a.bonds.count == 3 and b.b.bonds.count == 3:
            b.order = 4

    c1 = s.particles[1]
    c5 = s.particles[5]
    c1.linker = 'head'
    c5.linker = 'tail'

    if not is_capped:
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

    s.apply_forcefield(ff, charges='gasteiger')
    s.set_box(padding=10)

    sim = lmps.Simulation(s, print_to_screen=False, log='ps_mon_relax.log')
    sim.add(lmps.Init())
    sim.add_min(min_style='fire')
    sim.run()

    s.add_particle_bonding()
    return s


def polymer_chain(length, ff):
    return random_walk(monomer(ff), length, forcefield=ff)
