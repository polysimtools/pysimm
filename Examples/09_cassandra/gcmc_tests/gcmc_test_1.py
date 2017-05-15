from pysimm import system, cassandra

# Create blank system object (It can be formed in many ways here we basically need
#   a box with a certain size to populate it afterwards with gas molecules
sst = system.System()
bxSize = 100
sst.dim.dx = bxSize
sst.dim.dy = bxSize
sst.dim.dz = bxSize

# In order to run CASSANDRA GCMC one need to create the CASSANDRA object
css = cassandra.Cassandra(sst)

# Read the some default CASSANDRA .inp parameters file
gcmcPropsRead = css.read_input('default.inp')

# Set the gas to be purged in a box
gcmcPropsRead['adsorbers'] = [['ch4', 2000]]

gcmcPropsRead['Chemical_Potential_Info'] = -33.640934

# The prefix for the all files that will be created by this run
gcmcPropsRead['Run_Name'] = 'gcmc_test_1'
css.add_gcmc(**gcmcPropsRead)
css.run()

# css.system.visualize()
