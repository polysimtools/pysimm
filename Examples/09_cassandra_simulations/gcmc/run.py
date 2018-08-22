from pysimm import system, cassandra, lmps


def run(test=False):
    # In order to run CASSANDRA GCMC one need to create the CASSANDRA object
    sst = system.System()
    sst.dim = system.Dimension(dx=40, dy=40, dz=45, center=[0, 0, 0])
    sst.forcefield = 'trappe/amber'
    lmps.check_lmps_attr(sst)

    css = cassandra.Cassandra(sst)
    
    # Read the CASSANDRA .inp parameters file -- common way to setup simulations.
    # Any of the read properties can be modified here afterwards
    my_gcmc_props = css.read_input('props.inp')
    
    # The prefix for the all files that will be created by this run
    my_gcmc_props['Run_Name'] = 'gas_adsorb'
    
    # Set the gas (gas system) to be purged in a box
    specie1 = system.read_lammps('co2.lmps')
    specie2 = system.read_lammps('ch4.lmps')
    specie3 = system.read_lammps('m-xylene.lmps')
    for s in [specie1, specie2, specie3]:
        s.forcefield = 'trappe/amber'

    css.add_gcmc(species=[specie1, specie2, specie3],
                 is_rigid=[True, False, True],
                 max_ins=[2000, 1000, 500],
                 chem_pot=[-27.34, -29.34, -24.59],
                 out_folder='gas_adsorb_results', **my_gcmc_props)
    css.run()

    for pt in css.system.particle_types:
        pt.elem = pt.real_elem
    
    css.system.write_lammps('gas_adsorb.lmps')
    css.system.write_xyz('gas_adsorb.xyz')

if __name__ == '__main__':
    run(False)
