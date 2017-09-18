Example 6: Creating systems of multiple polymethyl methacrylate chains using the random walk application
========================================================================================================
by Dylan Anstine and Michael E. Fortunato

### Importing pysimm modules/packages

This example code creates a system of four polymethyl methacrylate chains using the random walk application. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield, the random_walk application from the pysimm.apps package, and the pmma monomer function from the pysimm models database.

```
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk
from pysimm.models.monomers.gaff2.pmma import monomer
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Creating a polymethyl methacrylate monomer

We'll create a monomer models for polymethyl methacrylate using the pysimm models database.

`pmma = monomer()`

### Creating a reference to a forcefield object

We will need to use a force field object more than once, so we create a reference to it by instantiating the object and storing a reference to it in variable **f**.

`f = forcefield.Dreiding()`

### Creating a system of polymer chains

The force field objects retrieved from the **forcefield.Dreiding** object contain both Lennard-Jones and Buckingham potential parameters. By default the pair style for non bonded interactions is set to buckingham, however the Lennard-Jones potential is a better choice during relaxation of new polymer bonds. We can switch between the two potentials by modifying the **pair_style** attribute of our monomer system **pmma**.

```
pmma.pair_style = 'lj'
```

The first system we will create will have four polymer chains of the same length. A quick and easy way to accomplish this is to build one polymer chain, and then replicate it 3 more times, instead of building 4 chains. To do this we will first use the random_walk function to build a polymer of 5 repeat units.

```
polymer = random_walk(pmma , nmon=5, forcefield=f)
```

Then we will use the **system.replicate()** function to create a system of 4 of these chains. The **rand** argument will perform a rigid random translation and rigid rotation within the simulation box, but does not change the polymer geometry. The replicate function doesn't check for overlaps when inserting molecules, so we use a low density here to hopefully avoid any overlaps. We'll save this new polymer, which has a uniform distribution of molecular weights, to the variable name uniform_polymer.

```
uniform_polymer = system.replicate(polymer, 4, density = 0.022, rand=True)
```

The next system of polymer chains we will build will contain 4 chains, each with a different number of repeat units. As a small example, we'll build chains of 2, 4, 6, and 8 repeat units respectively. The first chain we build will be done similarly to how we build our chain before. This simulation box this first chain is built in will be the same as the one we end up with with all 4 chains, so we use a density of 0.3/4 to end up with a final density around 0.3. We'll save the new polymer chain in the variable name nonuniform_polymer.

```
nonuniform_polymer = random_walk(pmma , nmon=2, forcefield=f, density=0.3/4)
```

Next, we will build three more chains into the same simulation box. To accomplish this we use the keyword argument **s_** which represents the system in which we will built our first chain.

```
nonuniform_polymer = random_walk(pmma , nmon=4, s_=nonuniform_polymer, forcefield=f)
nonuniform_polymer = random_walk(pmma , nmon=6, s_=nonuniform_polymer, forcefield=f)
nonuniform_polymer = random_walk(pmma , nmon=8, s_=nonuniform_polymer, forcefield=f)
```

The **system.System** objects have a class method **set_mm_dist()** to calculate molecular weight distributions in the system. If we calculate the molecular weight distributions of our uniform and nonuniform polymer systems, we should see a dispersity, or uniformity, of 1, and 1.2 respectively. The **set_mm_dist** function calculates the number average and weight average molecular weights, and the dispersity. We'll calculate these and then output the information to the screen.

```
uniform_polymer.set_mm_dist()
nonuniform_polymer.set_mm_dist()
```