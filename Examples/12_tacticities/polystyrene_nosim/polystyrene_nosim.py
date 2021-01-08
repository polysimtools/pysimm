from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk_tacticity
from pysimm.models.monomers.dreiding.ps import monomer as ps


def run(test=False):
    # we'll take a polystyrene monomer from the pysimm models database
    monomer = ps(is_capped=True)

    # we'll decorate atoms of PS monomer with tags that are needed for more elaborate polymerisation
    # procedure that takes into account tacticity
    setattr(monomer.particles[9], 'linker', 'mirror')
    setattr(monomer.particles[9], 'rnd_wlk_tag', 'head_cap')
    setattr(monomer.particles[13], 'rnd_wlk_tag', 'tail_cap')

    # we'll instantiate a Dreiding forcefield object for use later
    f = forcefield.Dreiding()
    
    # the monomers do not have any charges, so we will derive partial charges using the gasteiger algorithm
    monomer.apply_charges(f, charges='gasteiger')
    
    # the buckingham potential isn't great at small distances, and therefore we use the LJ potential while growing the polymer
    monomer.pair_style = 'lj/cut'
    
    # run the polymer random walk tacticity method with 10 total repeat units
    polymer = random_walk_tacticity(monomer, 10, forcefield=f, tacticity='isotactic', sim='no')

    # write a few different file formats
    polymer.write_xyz('isotactic.xyz')
    polymer.write_yaml('isotactic.yaml')
    polymer.write_lammps('isotactic.lmps')
    polymer.write_chemdoodle_json('isotactic.json')

    # run the polymer random walk tacticity method with 10 total repeat units
    polymer = random_walk_tacticity(monomer, 10, forcefield=f, tacticity='syndiotactic', sim='no')

    # write a few different file formats
    polymer.write_xyz('syndiotactic.xyz')
    polymer.write_yaml('syndiotactic.yaml')
    polymer.write_lammps('syndiotactic.lmps')
    polymer.write_chemdoodle_json('syndiotactic.json')
    
    # run the polymer random walk tacticity method with 10 total repeat units
    polymer = random_walk_tacticity(monomer, 10, forcefield=f, tacticity=0.2, sim='no')

    # write a few different file formats
    polymer.write_xyz('iso_20_syndio_80.xyz')
    polymer.write_yaml('iso_20_syndio_80.yaml')
    polymer.write_lammps('iso_20_syndio_80.lmps')
    polymer.write_chemdoodle_json('iso_20_syndio_80.json')


if __name__ == '__main__':
    run()
