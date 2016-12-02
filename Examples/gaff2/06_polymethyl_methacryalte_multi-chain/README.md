Example 6: Creating a system of multiple polymethyl methacrylate chains using the random walk application
=========================================================================================================

### Importing pysimm modules/packages

This example code creates a system of four polymethyl methacrylate chains using the random walk application. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield, and the random_walk application from the pysimm.apps package.

```
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Creating a polymethyl methacrylate monomer

[PubChem](https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds) offers a database of compounds accessible through a RESTful API. pysimm utilizes this API and allows users to create **system.System** objects from a puchem SMILES query. In this example, we will use the SMILES string "cc(C)(C(=O)OC)" to generate a polymethyl methacrylate monomer using the **system.read_puchem_smiles()** method. The smiles format supports radicals, which we will exploit here so that each carbon only has two hydrogens and can participate in new polymer bonds. The function makes an http request to the PubChem server, which returns a mol file. This mol file is interpreted and the function returns a **system.System** object that we store in variable **pmma**. This system now contains elemental composition bond connectivity, and bond orders.

`pmma = system.read_pubchem_smiles('cc(C)(C(=O)OC)')`

Additionally, we need to identify the two carbon atoms as linker atoms for polymerization. In this case, the atoms can be accessed by their tag. Visualization software can be very helpful if the tags of certain atoms are not known ahead of time. In this case we know that the carbon atoms those with tags 3 and 6, and we assign 'head' and 'tail' to their **linker** attributes.

```
pmma.particles[3].linker='head'
pmma.particles[6].linker='tail'
```

### Creating a reference to a forcefield object

We will need to use a force field object more than once, so we create a reference to it by instantiating the object and storing a reference to it in variable **f**.

`f = forcefield.Gaff2()`

### Automatically applying force field parameters

Our **system.System** object **pmma** contains a class method **apply_forcefield()** that accepts a **forcefield.Forcefield** object and automatically assigns force field parameters to our systems. In this example we use our previously created **forcefield.Gaff2** object and pass it to the **apply_forcefield()** function, as well as specify we would like to derive Gasteiger charges.

```
pmma.apply_forcefield(f, charges='gasteiger')
```

### Optimizing the monomer structure using LAMMPS

We'll use the fire algorithm to perofrm energy minimization of our systems. The **lmps** module in pysimm contains convenience methods to configure and execute a simulation. In this case we will use **lmps.quick_min()** passing our system object **pmma** and the min_style we want to use.

```
lmps.quick_min(pmma, min_style='fire')
```

### Creating a system of polymer chains

The force field objects retrieved from the **forcefield.Gaff2** object contain both Lennard-Jones and Buckingham potential parameters. By default the pair style for non bonded interactions is set to buckingham, however the Lennard-Jones potential is a better choice during relaxation of new polymer bonds. We can switch between the two potentials by modifying the **pair_style** attribute of our monomer system **pmma**.

```
pmma.pair_style = 'lj'
```

The **random_walk** application requires a reference monomer, the number of repeat units in the desired chain, and the force field from which new force field types will be acquired. By default, a new **system.System** object will be generated where the polymer chain will grow, and this system will have a density of 0.3 g/ml. We are going to grow four chains in this system, so we explicitly set the density to 0.3/4 g/ml. During polymerization, new force field terms will be determined automatically, and LAMMPS simulations will be performed to relax new polymer bonds. The random_walk function returns the newly created polymer **system.System** object.

`polymer = random_walk(pmma , nmon=5, forcefield=f, density=0.3/4)`

Now we build three more polymer chains, and indicate we want these chains to be build in our already existing polymer system **polymer**, by using the keyword argument **s_**.

```
polymer = random_walk(pmma , nmon=5, s_=polymer, forcefield=f)
polymer = random_walk(pmma , nmon=5, s_=polymer, forcefield=f)
polymer = random_walk(pmma , nmon=5, s_=polymer, forcefield=f)
```

The **copolymer** method requires a list of reference monomers, the total number of repeat units in the desired chain, and the force field from which new force field types will be acquired. By default the pattern the monomers added to the chain will iterate through the list of reference monomers. This is equivalent to setting the pattern to a list of 1s with the same length as the list of reference monomers. During polymerization, new force field terms will be determined automatically, and LAMMPS simulations will be performed to relax new polymer bonds. To change the number of processors used during simulation, the the **settings** keyword can be used to provide a dictionary of simulation settings. The random_walk function returns the newly created polymer **system.System** object.

`polymer = copolymer([pe, ps], 10, pattern=[1, 1], forcefield=f, settings={'np': 2})`

### Writing the polymer system to various file formats

The **system.System** class has various method to format our system data into different file types. Here we write xyz, yaml, lammps data, and chemdoodle json file formats.

```
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')
```