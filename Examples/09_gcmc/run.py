import os
from pysimm import system, cassandra

def run(test=False):
    # In order to run CASSANDRA GCMC one need to create the CASSANDRA object
    sst = system.System()
    sst.dim = system.Dimension(dx=45, dy=45, dz=45, center=[0, 0, 0])
    sst.forcefield = 'trappe/amber'
    css = cassandra.Cassandra(sst)
    
    # Read the CASSANDRA .inp parameters file -- common way to setup simulations.
    # Any of the read properties can be modified here afterwards
    my_gcmc_props = css.read_input(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'props.inp'))
    
    # The prefix for the all files that will be created by this run
    my_gcmc_props['Run_Name'] = 'gas_adsorb'
    
    # Set the gas (gas system) to be purged in a box
    specie1 = system.read_lammps(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'co2.lmps'))
    specie2 = system.read_lammps(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ch4.lmps'))
    specie3 = system.read_lammps(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'm-xylene.lmps'))
    
    css.add_gcmc(species=[specie1, specie2, specie3],
                 max_ins=[2000, 1000, 500],
                 chem_pot=[-27.34, -29.34, -24.59],
                 out_folder='gas_adsorb_results', **my_gcmc_props)
    css.run()
    
    
    for pt in css.final_sst.particle_types:
        pt.elem = pt.real_elem
    
    css.final_sst.write_lammps('gas_adsorb.lmps')
    css.final_sst.write_xyz('gas_adsorb.xyz')

if __name__ == '__main__':
    run(False)