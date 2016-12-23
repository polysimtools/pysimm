Example 3: Creating benzene and fixing incorrect bond orders
============================================================
by Michael E. Fortunato

### Importing pysimm modules/packages

This example code creates a benzene molecule using the PubChem compound search API. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield.

`from pysimm import system, lmps, forcefield`

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Searching PubChem for a compound

[PubChem](https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds) offers a database of compounds accessible through a RESTful API. pysimm utilizes this API and allows users to create **system.System** objects from a puchem SMILES or CID query. In this example, we will use the SMILES string "c1=cc=cc=c1" to generate a system using the **system.read_puchem_smiles()** method. The function makes an http request to the PubChem server, which returns a mol file. This mol file is interpreted and the function returns a **system.System** object that we store in variable **s**. This system now contains elemental composition bond connectivity, and bond orders.

`s = system.read_pubchem_smiles('c1=cc=cc=c1')`

### Modifying bond orders

Unfortunately, the bond orders in the mol file returned by the PubChem API represent alternating double bonds. When we assign force field parameters, we want the bonds between carbon atoms in the ring to be considered as aromatic bonds. To fix this, we can iterate through the bonds in our **system.System** object **s** by using the iterable **s.bonds** container, check if the two particles involved in the bond **a** and **b** are carbon atoms by accessing their **elem** attribute, and if so, change the **order** attribute of the bond to "A".

```
for b in s.bonds:
    if b.a.elem=='C' and b.b.elem=='C':
        b.order='A'
```

### Automatically applying force field parameters

Our **system.System** object **s** contains a class method **s.apply_forcefield()** that accepts a **forcefield.Forcefield** object and automatically assigns force field parameters to our system. In this example we instantiate a **forcefield.Gaff2** object and pass it to the **s.apply_forcefield()** function.

`s.apply_forcefield(forcefield.Gaff2())`

### Optimizing the structure using LAMMPS

We'll use a two stage minimization procedure, first using a steepest decent algorithm followed by a conjugate gradient algorithm. The **lmps** module in pysimm contains convenience methods to configure and execute a simulation. In this case we will use **lmps.quick_min()** passing our system object **s**, the min_style we want to use, and give each simulation a name so that the log files have identifiable names.

```
lmps.quick_min(s, min_style='sd', name='min_sd')
lmps.quick_min(s, min_style='cg', name='min_cg')
```

### Writing out a system to various file formats

The **system.System** class has various method to format our system data into different file types. Here we write xyz, yaml, lammps data, and chemdoodle json file formats.

```
s.write_xyz('benzene.xyz')
s.write_yaml('benzene.yaml')
s.write_lammps('benzene.lmps')
s.write_chemdoodle_json('benzene.json')
```