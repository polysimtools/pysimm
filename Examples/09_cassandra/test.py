from pysimm import system, lmps, cassandra

xyzFile = '_testInput/result.xyz'
mcfFile_poly = '_testInput/result.mcf'


# Test of the .inp file write
gcmcProps = {'Run_Name': '_testOutput/result.out',
             'Chemical_Potential_Info': -35.6136,
             'Molecule_Files': {'file1': [mcfFile_poly, 1],
                                'file2': ['_testInput/ch4.mcf', 1000]},
             'Fragment_Files': {'file1': ['_testInput/ch4.dat', 1]},
             'Start_Type': {'start_type':'read_config',
                            'species':[1, 0],
                            'file_name': xyzFile},
             'Property_Info': {'prop1': 'Density',
                               'prop2': 'Energy_Total',
                               'prop3': 'Nmols',
                               'prop4': 'Volume',
                               'prop5': 'Pressure'},
             'Box_Info': {'box_size': 175.112}
             }

gcmc = cassandra.GCMC(**gcmcProps)

gcmc.write('result.inp')

sst = lmps.read_lammps('Tilanga/data.lmps')

cs = cassandra.Cassandra(sst)
#sst.write_xyz(xyzFile)

#cs.readParams('Tilanga/gcmc.inp')
#cs.write_chk('_testOutput/result.chk')

cs.write_mcf(mcfFile_poly)


#css = cs.simulation()
#css.run()