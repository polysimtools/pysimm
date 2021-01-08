import sys

from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk


def monomer(**kwargs):
    isomer = kwargs.get('isomer', {})
    try:
        import os
        if isomer == 'endo':
            s = system.read_mol(
                os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'NbTMS_H2_endo.mol'))
        elif isomer == 'exo_endo':
            s = system.read_mol(
                os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'NbTMS_H2_exo_endo.mol'))
        else:
            s = system.read_mol(
                os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'NbTMS_H2_exo_exo.mol'))
    except:
        print("!!!could not read mol!!!")
        raise

    f = forcefield.Dreiding()

    for b in s.bonds:
        if b.a.bonds.count == 3 and b.b.bonds.count == 3:
            b.order = 4

    s.apply_forcefield(f)

    h = s.particles[7]
    t = s.particles[8]
    m = s.particles[32]

    s.remove_spare_bonding()

    s.set_box(padding=10)

    h.linker = 'head'
    t.linker = 'tail'
    m.linker = 'mirror'
    #	lmps.quick_min(s, min_style='fire')

    s.add_particle_bonding()
    return s


def polymer_chain(length):
    mon = monomer()
    polym = random_walk(mon, length, forcefield=forcefield.Dreiding())
    return polym
