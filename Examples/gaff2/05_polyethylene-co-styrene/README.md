Example 5: Creating a polyethylene-co-styrene chain using the random walk application
=====================================================================================

### Importing pysimm modules/packages

This example code creates a polyethylene-co-styrene chain using the random walk application. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield, and the copolymer application from the pysimm.apps.random_walk module.

```
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import copolymer
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Creating a polyethylene monomer

[PubChem](https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds) offers a database of compounds accessible through a RESTful API. pysimm utilizes this API and allows users to create **system.System** objects from a puchem SMILES query. In this example, we will use the SMILES string "cc" to generate a polyethylene monomer using the **system.read_puchem_smiles()** method. The smiles format supports radicals, which we will exploit here so that each carbon only has two hydrogens and can participate in new polymer bonds. The function makes an http request to the PubChem server, which returns a mol file. This mol file is interpreted and the function returns a **system.System** object that we store in variable **pe**. This system now contains elemental composition bond connectivity, and bond orders.

`pe = system.read_pubchem_smiles('cc')`

Additionally, we need to identify the two carbon atoms as linker atoms for polymerization. In this case, the atoms can be accessed by their tag. Visualization software can be very helpful if the tags of certain atoms are not known ahead of time. In this case we know that the carbon atoms are the first two atoms in our system, and we assign 'head' and 'tail' to their **linker** attributes.

```
pe.particles[1].linker='head'
pe.particles[2].linker='tail'
```

### Creating a polystyrene monomer

This time we use the SMILES string "cc(C1=CC=CC=C1)" to generate a polystyrene monomer using the **system.read_puchem_smiles()** method, and store the new system object in variable **ps**.

`ps = system.read_pubchem_smiles('cc')`

The linker atoms in this monomer are 7 and 8.

```
ps.particles[8].linker='head'
ps.particles[7].linker='tail'
```

Like in Example 3 where we create a benzene molecule, we need to fix the bond orders stored in our system object to reflect the aromaticity in the ring. This time we can't simply check to see if the particles in the bonds are carbons, we additionally need to ensure they are not linker atoms.

```
for b in ps.bonds:
    if (not b.a.linker and not b.b.linker) and b.a.elem=='C' and b.b.elem=='C':
        b.order='A'
```

### Creating a reference to a forcefield object

We will need to use a force field object more than once, so we create a reference to it by instantiating the object and storing a reference to it in variable **f**.

`f = forcefield.Gaff2()`

### Automatically applying force field parameters

Our **system.System** objects **pe** and **ps** contain class methods **s.apply_forcefield()** that accepts a **forcefield.Forcefield** object and automatically assigns force field parameters to our systems. In this example we use our previously created **forcefield.Gaff2** object and pass it to the **apply_forcefield()** function, as well as specify we would like to derive Gasteiger charges.

```
pe.apply_forcefield(f, charges='gasteiger')
ps.apply_forcefield(f, charges='gasteiger')
```

### Optimizing the monomer structure using LAMMPS

We'll use the fire algorithm to perofrm energy minimization of our systems. The **lmps** module in pysimm contains convenience methods to configure and execute a simulation. In this case we will use **lmps.quick_min()** passing our system object **pe** or **ps** and the min_style we want to use.

```
lmps.quick_min(pe, min_style='fire')
lmps.quick_min(ps, min_style='fire')
```

### Creating a copolymer chain

The force field objects retrieved from the **forcefield.Gaff2** object contain both Lennard-Jones and Buckingham potential parameters. By default the pair style for non bonded interactions is set to buckingham, however the Lennard-Jones potential is a better choice during relaxation of new polymer bonds. We can switch between the two potentials by modifying the **pair_style** attribute of our monomer systems **pe** and **ps**.

```
pe.pair_style = 'lj'
ps.pair_style = 'lj'
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