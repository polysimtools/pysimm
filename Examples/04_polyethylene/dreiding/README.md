Example 4: Creating a polyethylene chain using the random walk application
==========================================================================
by Michael E. Fortunato and Dylan Anstine

### Importing pysimm modules/packages

This example code creates a polyethylene chain using the random walk application. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield, the random_walk application from the pysimm.apps package, and the polyethylene monomer from the pysimm models database.

```
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk
from pysimm.models.monomers.dreiding.pe import monomer
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Creating a polyethylene monomer

The pysimm models database has a few convenience functions to build monomer repeat units that are ready to be polymerized. In this example we'll create a polyethylene monomer from the database. The monomer already contained information about which particles are linkers for polymerization, however there are no partial charges assigned to particles in these monomer models. We'll derive charges in the next step.

`pe = monomer()`

### Creating a reference to a forcefield object

We will need to use a force field object more than once, so we create a reference to it by instantiating the object and storing a reference to it in variable **f**.

`f = forcefield.Dreiding()`

### Deriving Gasteiger partial charges

Our **system.System** object **pe** contains a class method **pe.apply_charges()** that accepts a **forcefield.Forcefield** object and the name of the charge derivation method we would like to use.

`pe.apply_forcefield(f, charges='gasteiger')`

### Creating a polymer chain

The **random_walk** application requires a reference monomer, the number of repeat units in the desired chain, and the force field from which new force field types will be acquired. During polymerization, new force field terms will be determined automatically, and LAMMPS simulations will be performed to relax new polymer bonds. The random_walk function returns the newly created polymer **system.System** object. Note: The random walk application replicates monomers as they are in the reference system, meaning the ends of the resulting polymer chains have uncapped, dangling bonds.

`polymer = random_walk(pe, 10, forcefield=f)`

### Writing the polymer system to various file formats

Like above, we can write our system to various file formats.

```
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')
```

### Restarting a polymerization

If we want to pick up where we might have left off with a previous polymerization, we can restart by using the **copolymer** function in the random_walk application.  The random_walk application (including copolymer) uses the last head and tail particles in your growing chain. We will exploit this by calling our original polymer chain the first "monomer", and grow on the end of the chain.

First we need to import the copolymer function from the random_walk application.

`from pysimm.apps.random_walk import copolymer`

Now we can read in the yaml file we saved after building our first polymer chain of 10 units.

`original_polymer = system.read_yaml('polymer.yaml')`

The yaml file format tries to represent the system object as closely as possible, and retain the linker information the random_walk application needs. Now we can call the copolymer function, and give it a list of reference "monomers", although the first reference we pass will be our original polymer. The next argument is the *total* number of monomers we want to add. Using a value of 6 will ultimately mean we will insert 1 original_polymer and then 5 new monomers. This pattern is defined using the **pattern** argument. A list of [1, 5] means the first reference will be inserted once, the second reference will be inserted 5 times, and this pattern will be repeated until we reach the total number of insertions we defined in the second positional argument. Like before, we provide a forcefield object to acquire new forcefield types, although at this point we have everything we need on our system object.

`longer_polymer = copolymer([original_polymer, pe], 6, pattern=[1, 5], forcefield=f)`