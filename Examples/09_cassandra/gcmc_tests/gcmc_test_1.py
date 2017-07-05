from pysimm import system, cassandra
from collections import OrderedDict

bxSize = 44.2325

# In order to run CASSANDRA GCMC one need to create the CASSANDRA object
css = cassandra.Cassandra()

# Read the some default CASSANDRA .inp parameters file
gcmcPropsRead = css.read_input('my_props.inp')
gcmcPropsRead['Box_Info'] = OrderedDict([('box_count', 1), ('box_size', bxSize)])

# The prefix for the all files that will be created by this run
gcmcPropsRead['Run_Name'] = 'no_pim_box'
gcmcPropsRead['Rcutoff_Low'] = 0.

# Set the gas (gas system) to be purged in a box
h2 = system.read_lammps('h2.lmps')
gas_syst = cassandra.McSystem(h2, max_ins=5000, chem_pot=-23.34)

job = cassandra.GCMC(gas_syst, out_folder='gcmc_test_1', **gcmcPropsRead)
css.add_gcmc(job)

css.run()

job.tot_sst.write_lammps('hydrogens.lmps')
job.tot_sst.system.visualize()
