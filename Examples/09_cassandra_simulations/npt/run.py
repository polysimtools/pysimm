from pysimm import system, cassandra, lmps


def run(test=False):
    # Setup the empty cubic box with
    bx_size = 60
    sst = system.System()
    sst.dim = system.Dimension(dx=bx_size, dy=bx_size, dz=bx_size, center=[bx_size / 2, bx_size / 2, bx_size / 2])
    sst.forcefield = 'trappe/amber'

    molec = system.read_lammps('c2h4.lmps')
    molec.forcefield = 'trappe/amber'

    cs = cassandra.Cassandra(sst)
    npt_props = cs.read_input('props.inp')
    npt_props['Pressure_Info'] = 25  # Simulated pressure in bars
    npt_props['Start_Type'] = {'start_type': 'make_config', 'species': 300}
    npt_props['Simulation_Length_Info'] = {'run': 10000}
    npt_props['Property_Info'] = {'prop1': 'energy_total',
                                  'prop2': 'volume',
                                  'prop3': 'mass_density'}

    cs.add_simulation('NPT', species=molec, is_rigid=True, out_folder='results', **npt_props)

    cs.run()

    lmps.check_lmps_attr(cs.system)
    cs.system.write_lammps('final_conf.lmps')

if __name__ == '__main__':
    run(False)
