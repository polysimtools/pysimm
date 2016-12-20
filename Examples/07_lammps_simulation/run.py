from pysimm import system, lmps, forcefield

logFileName = 'steps.log'

# use an import from .yaml to get of previously created pmma-based 4-chain polymer structure
polymer = system.read_yaml("polymer.yaml")

# Initialize the wrapper object around the polymer that will organize the work with LAMMPS
print('Creating Simulation object with log-file in "{0:s}"\n'.format(logFileName))
sim = lmps.Simulation(polymer, log= 'steps.log')

# setting up the parameters for the energy optimization
#  add_min() method will add the "Minimization" task to the task que of the 
#  Simulation object that is stored in sim.sim list
print('Creating and running energy-minimization task:')
sim.add_min(min_style = 'fire', name = 'min_fire',  etol = 1.0e-5, ftol = 1.0e-5)

print('Creating and running molecular dynamics task:')
# Let's set up the Molecular dynamics task
sim.add_md(ensemble='nvt', timestep=0.5)

print('List of simulation tasks ready to run:')
print(sim.sim)

print('Input that will be passed to LAMMPS when the simulation is performed:')
sim.write_input()
print(sim.input)

# call run to run the simulation
sim.run()
