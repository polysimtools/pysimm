from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk_tacticity
from pysimm.models.monomers.dreiding.ps_tacticity import monomer as PS_tacticity

def run(test=False):
    # we'll make a polystyrene monomer from the pysimm models database
    A = PS_tacticity()
    
    # we'll instantiate a Dreiding forcefield object for use later
    f = forcefield.Dreiding()
    
    # the monomers do not have any charges, so we will derive partial charges using the gasteiger algorithm
    A.apply_charges(f, charges='gasteiger')
    
    # the buckingham potential isn't great at small distances, and therefore we use the LJ potential while growing the polymer
    A.pair_style = 'lj/cut'
    
    # run the polymer random walk tacticity method with 10 total repeat units
    polymer = random_walk_tacticity(A, 10, forcefield=f,capped=True,tacticity='isotactic',sim='no')

    # write a few different file formats
    polymer.write_xyz('isotactic.xyz')
    polymer.write_yaml('isotactic.yaml')
    polymer.write_lammps('isotactic.lmps')
    polymer.write_chemdoodle_json('isotactic.json')

    # run the polymer random walk tacticity method with 10 total repeat units
    polymer = random_walk_tacticity(A, 10, forcefield=f,capped=True,tacticity='syndiotactic',sim='no')

    # write a few different file formats
    polymer.write_xyz('syndiotactic.xyz')
    polymer.write_yaml('syndiotactic.yaml')
    polymer.write_lammps('syndiotactic.lmps')
    polymer.write_chemdoodle_json('syndiotactic.json')
    
    # run the polymer random walk tacticity method with 10 total repeat units
    polymer = random_walk_tacticity(A, 10, forcefield=f,capped=True,tacticity=0.2,sim='no')

    # write a few different file formats
    polymer.write_xyz('iso_20_syndio_80.xyz')
    polymer.write_yaml('iso_20_syndio_80.yaml')
    polymer.write_lammps('iso_20_syndio_80.lmps')
    polymer.write_chemdoodle_json('iso_20_syndio_80.json')
    
if __name__ == '__main__':
    run()
