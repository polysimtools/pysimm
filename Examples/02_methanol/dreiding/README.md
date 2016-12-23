Example 2: Creating methanol using the PubChem compound search API
=================================================================
by Michael E. Fortunato

### Importing pysimm modules/packages

This example code creates a methanol molecule using the PubChem compound search API. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield.

`from pysimm import system, lmps, forcefield`

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Searching PubChem for a compound

[PubChem](https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds) offers a database of compounds accessible through a RESTful API. pysimm utilizes this API and allows users to create **system.System** objects from a puchem SMILES or CID query. In this example, we will use the SMILES string "CO" to generate a system using the **system.read_puchem_smiles()** method. The function makes an http request to the PubChem server, which returns a mol file. This mol file is interpreted and the function returns a **system.System** object that we store in variable **s**. This system now contains elemental composition bond connectivity, and bond orders.

`s = system.read_pubchem_smiles('CO')`

### Automatically applying force field parameters

Our **system.System** object **s** contains a class method **s.apply_forcefield()** that accepts a **forcefield.Forcefield** object and automatically assigns force field parameters to our system. In this example we instantiate a **forcefield.Dreiding** object and pass it to the **s.apply_forcefield()** function.

`s.apply_forcefield(forcefield.Dreiding())`

### Optimizing the structure using LAMMPS

We'll use the fire algorithm to perofrm energy minimization of our system. The **lmps** module in pysimm contains convenience methods to configure and execute a simulation. In this case we will use **lmps.quick_min()** passing our system object **s** and the min_style we want to use.

`lmps.quick_min(s, min_style='fire')`

### Writing out a system to various file formats

The **system.System** class has various method to format our system data into different file types. Here we write xyz, yaml, lammps data, and chemdoodle json file formats.

```
s.write_xyz('methanol.xyz')
s.write_yaml('methanol.yaml')
s.write_lammps('methanol.lmps')
s.write_chemdoodle_json('methanol.json')
```