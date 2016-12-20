Example X: Creating a system of a mixture of ethanol and acetone, and run NPT simulation using LAMMPS
=========================================================================================================

### Importing pysimm modules/packages

This example code creates a system of mixture of ethanol and acetone using **system.replicate** command. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield.

```
from pysimm import system, lmps, forcefield
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Creating a ethanol and acetone molecule

[PubChem](https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds) offers a database of compounds accessible through a RESTful API. pysimm utilizes this API and allows users to create **system.System** objects from a puchem SMILES query. In this example, we will use the SMILES string "CCO" to generate a ethanol molecule using the **system.read_puchem_smiles()** method. The function makes an http request to the PubChem server, which returns a mol file. This mol file is interpreted and the function returns a **system.System** object that we store in variable **pmma**. This system now contains elemental composition bond connectivity, and bond orders.

`ethanol = system.read_pubchem_smiles('CCO')`

Next, we will create a acetone molecule object the same way.

`acetone = system.read_pubchem_smiles('CC(=O)C')`

### Creating a reference to a forcefield object

We will need to use a force field object more than once, so we create a reference to it by instantiating the object and storing a reference to it in variable **f**.

`f = forcefield.Gaff2()`

### Automatically applying force field parameters

Our **system.System** object **ethanol** and **acetone** contains a class method **apply_forcefield()** that accepts a **forcefield.Forcefield** object and automatically assigns force field parameters to our systems. In this example we use our previously created **forcefield.Gaff2** object and pass it to the **apply_forcefield()** function, as well as specify we would like to derive Gasteiger charges.

```
ethanol.apply_forcefield(f, charges='gasteiger')
acetone.apply_forcefield(f, charges='gasteiger')
```
If **Ambertools** is installed, and **antechamber** is in the PATH, the charges can also assigned using AM1-BCC charge method, which is more compatible with GAFF2 force field. The corresponding commands are:

```
amber.calc_charges(ethanol)
amber.calc_charges(acetone)
```

### Optimizing the monomer structure using LAMMPS

We'll use the fire algorithm to perofrm energy minimization of our systems. The **lmps** module in pysimm contains convenience methods to configure and execute a simulation. In this case we will use **lmps.quick_min()** passing our system object **ethanol**/**acetone** and the min_style we want to use.

```
lmps.quick_min(ethanol, min_style='fire')
lmps.quick_min(acetone, min_style='fire')
```

### Creating a system of a mixture of 300 ethanol and 200 acetone

The command to generate a system of 300 ethanol and 200 acetone is **system.replicate**, which requires a list containing the molecules and a list specifying the numbers of each molecule. 

```
molecule_list=[ethanol,acetone]
n_molecules=[300,200]
```

The **system.replicate** command will also take the density input to set up the simulation box. 

`s=system.replicate(molecule_list, n_molecules , density=0.3)`

Now that we have the molecules packed in a box, we can carry out simple minimization, NVT molecular dynamics for 10ps, and NPT molecular dynamics for 100ps to obtain a reasonable packed mixture of ethanol and acetone. 

```
lmps.quick_min(s, name='fire_min', min_style='fire', print_to_screen=True, thermo=500)
lmps.quick_md(s, name='nvt_md', ensemble='nvt', length=10000, print_to_screen=True, thermo=500)
lmps.quick_md(s, name='npt_md', ensemble='npt', length=100000, print_to_screen=True, thermo=500)
```

### Writing the mixture system to various file formats

The **system.System** class has various method to format our system data into different file types. Here we write xyz, yaml, lammps data, and chemdoodle json file formats.

```
s.write_xyz('mixture.xyz')
s.write_yaml('mixture.yaml')
s.write_lammps('mixture.lmps')
s.write_chemdoodle_json('mixture.json')
```