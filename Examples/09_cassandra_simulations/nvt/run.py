from pysimm import system, cassandra, lmps


def run(test=False):
    # Setup the box with acetelene molecules on the regular grid
    sst = system.System()

    bx_size = 30
    sst.dim = system.Dimension(dx=bx_size, dy=bx_size, dz=bx_size, center=[bx_size / 2, bx_size / 2, bx_size / 2])
    sst.forcefield = 'trappe/amber'

    molec = system.read_lammps('c2h4.lmps')
    molec.forcefield = 'trappe/amber'

    cs = cassandra.Cassandra(sst)
    nvt_props = cs.read_input('props.inp')

    nvt_props['Temperature_Info'] = 400
    nvt_props['Start_Type'] = {'start_type': 'make_config', 'species': 300}
    nvt_props['Simulation_Length_Info'] = {'run': 300000}
    nvt_props['Property_Info'] = {'prop1': 'energy_total',
                                  'prop2': 'pressure',
                                  'prop3': 'mass_density'}
    cs.add_nvt(species=molec, is_rigid=True, out_folder='results', **nvt_props)

    cs.run()

    lmps.check_lmps_attr(cs.system)
    cs.system.write_lammps('final_conf.lmps')

if __name__ == '__main__':
    run(False)
