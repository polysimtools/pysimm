from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk_tacticity
from pysimm.models.monomers.dreiding.NbTMS_H2_tacticity import monomer as NbTMS
from pysimm.apps.equilibrate import equil


def run(test=False):
    # we'll make a vinyl-type polynorbornene monomer from the pysimm models database
    A = NbTMS(isomer="exo_exo")

    setattr(A.particles[1], 'rnd_wlk_tag', 'tail_cap')
    setattr(A.particles[49], 'rnd_wlk_tag', 'head_cap')

    # we'll instantiate a Dreiding forcefield object for use later
    f = forcefield.Dreiding()

    # the monomers do not have any charges, so we will derive partial charges using the gasteiger algorithm
    A.apply_charges(f, charges='gasteiger')

    # the buckingham potential isn't great at small distances, and therefore we use the LJ potential while growing the polymer
    A.pair_style = 'lj/cut'

    # run the polymer random walk tacticity method with 10 total repeat units, all isotactic insertions, no minimizations
    polymer = random_walk_tacticity(A, 10, forcefield=f, capped=True, tacticity='isotactic', sim='no', rotation=180)
    write_file_formats(polymer, "isotactic_PNB_nosim10")

    # run the polymer random walk tacticity method with 10 total repeat units, all syndiotactic insertions, no minimizations
    polymer = random_walk_tacticity(A, 10, forcefield=f, capped=True, tacticity='syndiotactic', sim='no', rotation=180)
    write_file_formats(polymer, "syndiotactic_PNB_nosim10")

    # run the polymer random walk tacticity method with 50 total repeat units with minimizations and error
    # checking (when a hardcore overlap is found, the last inserted monomer shrunken and gradually re-expanded)
    polymer = random_walk_tacticity(A, 50, forcefield=f, capped=True, tacticity='syndiotactic', error_check=True,
                                    rotation=180, md_spacing=1)
    write_file_formats(polymer, "syndiotactic_PNB_sim50")

    # pack and equilibrate eight 20mers of exo_exo and exo_endo isomers to compare densities
    pack_and_equil(A, 20, 8, f, "exo")
    A = NbTMS(isomer="exo_endo")
    A.apply_charges(f, charges='gasteiger')
    A.pair_style = 'lj/cut'
    pack_and_equil(A, 20, 8, f, "endo")


def pack_and_equil(A, n, x, f, prefix):
    # A: monomer
    # n: number of monomers
    # x: number of chains
    # f: forcefield
    # prefix: prefix for output files

    # run the polymer random walk tacticity method with n total repeat units
    polymer = random_walk_tacticity(A, n, forcefield=f, capped=True, tacticity='syndiotactic', rotation=180,
                                    error_check=False, sim=0)
    write_file_formats(polymer, prefix + "_1", unwrap=True)

    # quick opt of polymer
    lmps.quick_min(polymer, min_style='cg', etol=1.0e-4, maxiter=100000)

    # write a few different file formats
    polymer.unwrap()
    write_file_formats(polymer, prefix + "_1_cg")

    # pack x copies of polymer
    polymers = system.replicate(polymer, x, density=0.005)
    # polymers = polymer
    write_file_formats(polymers, prefix + "_" + str(x))
    lmps.quick_min(polymers, min_style='cg', etol=1.0e-4, maxiter=100000)
    write_file_formats(polymers, prefix + "_" + str(x) + "_cg")

    # quickmd
    nvt_settings = {
        'name': 'nvt_md',
        'print_to_screen': True,
        'ensemble': 'nvt',
        'temperature': {
            'start': 100,
            'stop': 300
        },
        'new_v': True,
        'length': 10000
    }
    npt_settings = {
        'name': 'npt_md',
        'print_to_screen': True,
        'ensemble': 'npt',
        'temperature': 300,
        'new_v': True,
        'pressure': {
            'start': 1000,
            'stop': 1
        },
        'length': 100000,
        'thermo_style': 'custom step temp press density'
    }
    # npt calcs need "add neigh_modify" command to reneighbor more often during compression of npt step
    sim = lmps.Simulation(polymers, name='npt_reneighbor', debug=True)
    sim.add_custom('neigh_modify delay 0')
    sim.add(lmps.Velocity(temperature=1000))
    sim.add_md(length=10000, ensemble='npt', temperature=1000, pressure=5000)
    sim.run()
    write_file_formats(polymers, prefix + "_" + str(x) + "_npt")
    write_file_formats(polymers, prefix + "_" + str(x) + "_npt_unwrapped", unwrap=True)

    # 21-step equilibration
    equil(polymers, np=1, pmax=50000)
    write_file_formats(polymers, prefix + "_" + str(x) + "_equil")
    write_file_formats(polymers, prefix + "_" + str(x) + "_equil_unwrapped", unwrap=True)


def write_file_formats(s, name, **kwargs):
    unwrap = kwargs.get('unwrap', False)
    if unwrap:
        s.unwrap()
    s.write_xyz(name + '.xyz')
    s.write_yaml(name + '.yaml')
    s.write_lammps(name + '.lmps')
    s.write_chemdoodle_json(name + '.json')
    if unwrap:
        s.wrap()


if __name__ == '__main__':
    run()
