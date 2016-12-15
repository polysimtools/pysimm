Example 4: Creating a polyethylene chain using the random walk application
==========================================================================

### Importing pysimm modules/packages

This example code creates a polyethylene chain using the random walk application. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield, and the random_walk application from the pysimm.apps package.

```
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Creating a polyethylene monomer

[PubChem](https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds) offers a database of compounds accessible through a RESTful API. pysimm utilizes this API and allows users to create **system.System** objects from a puchem SMILES or CID query. In this example, we will use the SMILES string "cc" to generate a polyethylene monomer using the **system.read_puchem_smiles()** method. The smiles format supports radicals, which we will exploit here so that each carbon only has two hydrogens and can participate in new polymer bonds. The function makes an http request to the PubChem server, which returns a mol file. This mol file is interpreted and the function returns a **system.System** object that we store in variable **s**. This system now contains elemental composition bond connectivity, and bond orders.

`s = system.read_pubchem_smiles('cc')`

Additionally, we need to identify the two carbon atoms as linker atoms for polymerization. In this case, the atoms can be accessed by their tag. Visualization software can be very helpful if the tags of certain atoms are not known ahead of time. In this case we know that the carbon atoms are the first two atoms in our system, and we assign 'head' and 'tail' to their **linker** attributes.

```
s.particles[1].linker='head'
s.particles[2].linker='tail'
```

### Creating a reference to a forcefield object

We will need to use a force field object more than once, so we create a reference to it by instantiating the object and storing a reference to it in variable **f**.

`f = forcefield.Dreiding()`

### Automatically applying force field parameters

Our **system.System** object **s** contains a class method **s.apply_forcefield()** that accepts a **forcefield.Forcefield** object and automatically assigns force field parameters to our system. In this example we use our previously created **forcefield.Dreiding** object and pass it to the **s.apply_forcefield()** function, as well as specify we would like to derive Gasteiger charges.

`s.apply_forcefield(f, charges='gasteiger')`

### Optimizing the monomer structure using LAMMPS

We'll use the fire algorithm to perofrm energy minimization of our system. The **lmps** module in pysimm contains convenience methods to configure and execute a simulation. In this case we will use **lmps.quick_min()** passing our system object **s** and the min_style we want to use.

`lmps.quick_min(s, min_style='fire')`

### Writing the monomer system to various file formats

The **system.System** class has various method to format our system data into different file types. Here we write xyz, yaml, lammps data, and chemdoodle json file formats.

```
s.write_xyz('pe_monomer.xyz')
s.write_yaml('pe_monomer.yaml')
s.write_lammps('pe_monomer.lmps')
s.write_chemdoodle_json('pe_monomer.json')
```

### Creating a polymer chain

The **random_walk** application requires a reference monomer, the number of repeat units in the desired chain, and the force field from which new force field types will be acquired. During polymerization, new force field terms will be determined automatically, and LAMMPS simulations will be performed to relax new polymer bonds. The random_walk function returns the newly created polymer **system.System** object. Note: The random walk application replicates monomers as they are in the reference system, meaning the ends of the resulting polymer chains have uncapped, dangling bonds.

`polymer = random_walk(s, 10, forcefield=f)`

### Writing the polymer system to various file formats

Like above, we can write our system to various file formats.

```
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')
```