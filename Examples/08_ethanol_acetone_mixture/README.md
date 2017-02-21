Example 8: Creating a system of a mixture of ethanol and acetone, and run NPT simulation using LAMMPS
=========================================================================================================
by Ping Lin and Michael E. Fortunato

### Importing pysimm modules/packages

This example code creates a system representing a mixture of ethanol and acetone using **system.replicate** command. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield.

```
from pysimm import system, lmps, forcefield
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Creating an ethanol and acetone molecule

[PubChem](https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds) offers a database of compounds accessible through a RESTful API. pysimm utilizes this API and allows users to create **system.System** objects from a puchem SMILES query. In this example, we will use the SMILES string "CCO" to generate a ethanol molecule using the **system.read_puchem_smiles()** method. The function makes an http request to the PubChem server, which returns a mol file. This mol file is interpreted and the function returns a **system.System** object that we store in variable **ethanol**. This system now contains elemental composition bond connectivity, and bond orders.

`ethanol = system.read_pubchem_smiles('CCO')`

Next, we will create a acetone molecule object the same way.

`acetone = system.read_pubchem_smiles('CC(=O)C')`

### Creating a reference to a forcefield object

We will need to use a force field object more than once, so we create a reference to it by instantiating the object and storing a reference to it in variable **f**.

`f = forcefield.Gaff2()`

### Automatically applying force field parameters

Our **system.System** objects **ethanol** and **acetone** contain a class method **apply_forcefield()** that accepts a **forcefield.Forcefield** object and automatically assigns force field parameters to our systems. In this example we use our previously created **forcefield.Gaff2** object and pass it to the **apply_forcefield()** function, as well as specify we would like to derive Gasteiger charges.

```
ethanol.apply_forcefield(f, charges='gasteiger')
acetone.apply_forcefield(f, charges='gasteiger')
```
If **Ambertools** is installed, and **antechamber** is in the PATH, the charges can also assigned using AM1-BCC charge method, which is more compatible with GAFF2 force field. The corresponding commands are:

```
import amber
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

The **system.replicate** command will also take density input in g/ml to set up the simulation box. 

`s=system.replicate(molecule_list, n_molecules , density=0.3)`

Now that we have the molecules packed in a box, we can carry out simple minimization, NVT molecular dynamics for 10ps, and NPT molecular dynamics for 100ps to obtain a reasonable packed mixture of ethanol and acetone. The settings for each simulation step can be configured in a python dictionary.

```
min_settings = {
    'name': 'fire_min',
    'min_style': 'fire',
}
```

During the NVT simulation, we will generate new velocities and heat our system from 100 K to 300 K.

```
nvt_settings = {
    'name': 'nvt_md',
    'ensemble': 'nvt',
    't_start': 100,
    't_stop': 300,
    'new_v': True,
    'length': 10000
}
```

During the NPT simulation, we'll again generate new initial velocities and start with a positive compressive pressure and slowly decrease the pressure to 1 atm. We can also use the lammps keyword thermo_style to change the thermodynamics output in our log file. Here we'll output the step, temperature, pressure and density of our system.

```
npt_settings = {
    'name': 'npt_md',
    'ensemble': 'npt',
    'temp': 300,
    'new_v': True,
    'p_start': 1000,
    'p_stop': 1,
    'length': 100000,
    'thermo_style': 'custom step temp press density'
}
```

The simulations then can be performed using the **quick_min** and **quick_md** functions in the **lmps** module, and supplying the unpacked dictionary as keyword arguments.
```
lmps.quick_min(s, **min_settings)
lmps.quick_md(s, **nvt_settings)
lmps.quick_md(s, **npt_settings)
```

### Writing the mixture system to various file formats

The **system.System** class has various method to format our system data into different file types. Here we write xyz, yaml, lammps data, and chemdoodle json file formats.

```
s.write_xyz('mixture.xyz')
s.write_yaml('mixture.yaml')
s.write_lammps('mixture.lmps')
s.write_chemdoodle_json('mixture.json')
```