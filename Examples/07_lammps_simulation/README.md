Example 7: Running energy minimization and molecular dynamics simulations of a polymer system using LAMMPS
=======================================================================
by Alexander Demidov and Michael E. Fortunato

### Importing pysimm modules/packages
This example code performs a simulation using the LAMMPS software package performing both energy optimization and molecular dynamics  routines of 4 polymer chains created in Example 6. The structure file we need from the previous example is 'polymer.yaml'. The yaml file written from pysimm attempts to recreate the **system.System** data structure as closely as possible, and contains the force field information we will need when running our simulation. Alternatively, you could use the 'polymer.lmps' data file written from Example 6 as well, as this also contains all the information needed to begin our simulation. We've provided a sample 'polymer.yaml' file here, but feel free to use the polymer you created previously.

First we import the required modules.

```
from pysimm import system, lmps, forcefield
```

### Reading the yaml file containing system information
We will use the **system.read_yaml** function to read the yaml file we created in Example 6 in either the dreiding or gaff2 example section. Copy either the uniform_polymer or nonuniform_polymer yaml file here and name it polymer.yaml. We'll store a reference to the system in the variable name **polymer**.

```
polymer = system.read_yaml("polymer.yaml")
```

### Creating a Simulation object to define the simulation workflow

In order to perform either energy optimizations or molecular dynamics simulations we need to wrap the object of Simulation class around polymer system object **polymer**.  The Simulation object is needed to organize and structure the input to LAMMPS. LAMMPS will create a log file, and we can give it a name; here we call it 'steps.log'. We'll store a reference to our Simulation object in the variable **sim**.

```
sim = lmps.Simulation(polymer, log= 'steps.log')
```

### Adding tasks to the Simulation object

First we would like to perform an energy minimization of our polymer system. We use the **Simulation.add_min()** class method to add the task to our Simulation object **sim**. Ultimately we will add a **lmps.Minimization** to our Simulation object. The **add_min** method accepts either a **Minimization** object, or a set of key word arguments that will be passed to the **Minimization** constructor. Here we'll explicitly set the minimization style for the fire algorithm, the energy tolerance and force tolerance, etol and ftol, to 1e-5, and name our task min_fire. See documentation for a complete list of default values for **lmps.Minimization**.

```
sim.add_min(min_style = 'fire', name = 'min_fire',  etol = 1.0e-5, ftol = 1.0e-5)
```

Next let’s add to the queue of the sim object a molecular dynamics simulations task. Like before, ultimately we will be adding a **lmps.MolecularDynamics** object to our Simulation, however we can alternatively use key word arguments which will be passed to the MolecularDynamics constructor. Here we'll explicitly set the ensemble we would like to sample, and the timestep of our simulation. See documentation for a complete list of default values for **lmps.MolecularDynamics**.

```
sim.add_md(ensemble='nvt', timestep=0.5)
```

Now the sim object contains **2** tasks in the queue: first, the optimization and second is molecular dynamic. These are stored in the class attribute **sim**.

```
print(sim.sim)
```

To see the input that will be provided to LAMMPS, call the **Simulation.write_input()** function, and print the newly created attribute **Simulation.input**.

```
sim.write_input()
print(sim.input)
```

Finally, to run the simulation, call the **Simulation.run()** class method.

```
sim.run()
```
