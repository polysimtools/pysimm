Example 5: Creating a polyethylene-co-styrene chain using the random walk application
=====================================================================================
by Michael E. Fortunato

### Importing pysimm modules/packages

This example code creates a polyethylene-co-styrene chain using the random walk application. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield, the copolymer application from the pysimm.apps.random_walk module, and two monomer models from the pysimm models database.

```
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import copolymer
from pysimm.models.monomers.dreiding.pe import monomer as pe_monomer
from pysimm.models.monomers.dreiding.ps import monomer as ps_monomer
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Creating polyethylene and polystyrene monomers

We'll create two monomer models for polyethyene and polystyrene using the pysimm models database.

```
pe = pe_monomer()
ps = ps_monomer()
```

### Creating a reference to a forcefield object

We will need to use a force field object more than once, so we create a reference to it by instantiating the object and storing a reference to it in variable **f**.

`f = forcefield.Dreiding()`

### Deriving Gasteiger partial charges

**system.System** objects contains a class method **apply_charges()** that accepts a **forcefield.Forcefield** object and the name of the charge derivation method we would like to use. Here we'll derive gasteiger charges for both of our monomers.

```
pe.apply_charges(f, charges='gasteiger')
ps.apply_charges(f, charges='gasteiger')
```

### Creating a copolymer chain

The force field objects retrieved from the **forcefield.Dreiding** object contain both Lennard-Jones and Buckingham potential parameters. By default the pair style for non bonded interactions is set to buckingham, however the Lennard-Jones potential is a better choice during relaxation of new polymer bonds. We can switch between the two potentials by modifying the **pair_style** attribute of our monomer systems **pe** and **ps**. To switch back to the buckingham potential, reassign the **pair_style** attribute to 'buck'. This tells pysimm what parameters to write out for the simulation package later, as well as how to initialize the simulation.

```
pe.pair_style = 'lj'
ps.pair_style = 'lj'
```

The **copolymer** method requires a list of reference monomers, the total number of repeat units in the desired chain, and the force field from which new force field types will be acquired. For each reference monomer, the value in **pattern** at the same index will indicate how many monomers to add to the chain. For example, [1, 1] in this context will be an alternating copolymer. [1, 2] will be a pattern that has 1 pe monomer, followed by 2 ps monomers, repeated. [2, 5] will be a pattern that has 2 pe monomers, followed by 5 ps monomers, repeated, etc. During polymerization, new force field terms will be determined automatically, and LAMMPS simulations will be performed to relax new polymer bonds. The random_walk function returns the newly created polymer **system.System** object.

`polymer = copolymer([pe, ps], 10, pattern=[1, 1], forcefield=f)`

### Writing the polymer system to various file formats

The **system.System** class has various method to format our system data into different file types. Here we write xyz, yaml, lammps data, and chemdoodle json file formats.

```
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')
```