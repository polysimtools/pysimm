Example 9.1: The Grand-Canonical Monte Carlo simulations using Cassandra
========================================================================
by Alexander Demidov and Michael E. Fortunato

### Importing pysimm modules/packages

The example illustrates the basic principles of utilizing the **pysimm.cassandra** module. 
The code sets up the Grand-Canonical (&#956;VT) Monte Carlo simulations to be performed by the  [Cassandra](https://cassandra.nd.edu) software. 
The example shows how to use simultaneously molecules of different types and save the simulation results.

First, let's import the required pysimm modules: system and cassandra

```python
from pysimm import system, cassandra
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from GitHub is in your PYTHONPATH. 
See installation instructions for further directions.

### Creating the pysimm system

The cassandra module of the pysimm works with **pysimm.system** object modifying it (by addition or replacement) with results of Cassandra simulations. 
For this example the initial system will be an empty cubic cell with the size of 45 &#8491; centered at the origin of the coordinate system. 
Additionally, we explicitly specify the forcefield model for all particles that will be inserted into the system 

```python
sst = system.System()
sst.dim = system.Dimension(dx=45, dy=45, dz=45, center=[0, 0, 0])
sst.forcefield = 'trappe/amber'
```

The pysimm system is passed to the constructor of the main cassandra simulation object 

```python
css = cassandra.Cassandra(sst)
```

### Creating the simulation properties

The properties of the Cassandra simulation for the pysimm.cassandra module are provided in the form of a dictionary in which the key is the property name and the value (depending on the property itself) can be integer, float, string or another dictionary. 
Thus, all required Cassandra simulation properties can be specified directly in the Python run file by creating the dictionary. 
Alternatively, pysimm.cassandra object can read and interpret the proper Cassandra input file. 
If the files do not specify all properties required for &#956;VT simulations, they will be set to the default ones. 

Let's read the data from the **props.inp** file that is in the same example directory 

```python
my_gcmc_props = css.read_input('props.inp')
```
Now different simulation properties can be modified if needed. Let's, for example, change the prefix of all simulation files that will be created by Cassandra

```python
my_gcmc_props['Run_Name'] = 'gas_adsorb'
```

### Specifying the gas molecules

The setup of a Monte Carlo simulations should include also the information about entities to be inserted (modified or deleted) by the Monte Carlo algorithm. 
In our example, those are the gas molecules of three different species (CO<sub>2</sub>, methane, and m-xylene). 
For the **pysimm.cassandra** module they should be provided in the form of the **pysimm.system** objects.

The example folder contains the LAMMPS files that describe one molecule of each gas species, 
so the corresponding  **pysimm.system** objects can be created by reading the files by the **system.read_lammps()** method.

```python
specie1 = system.read_lammps('co2.lmps')
specie2 = system.read_lammps('ch4.lmps')
specie3 = system.read_lammps('m-xylene.lmps')
```

### Setup and run the simulations

The set-up of GCMC simulations in the **pysimm.cassandra** module is implemented in GCMC class. 
The way to add it to the simulation queue is to call **css.add_gcmc()** method. The method has few required keyword arguments 
* *species*: the molecules to be inserted;
* *chem_pot*: the values of chemical potential for each of species in the system;
* *max_ins*: maximal possible number of molecules of each species in the system;
* *out_folder*: relative path for all Cassandra files such as log files, .xyz, and .chk files

```python
css.add_gcmc(species=[specie1, specie2, specie3],
             max_ins=[2000, 1000, 500],
             chem_pot=[-27.34, -29.34, -24.59],
             out_folder='gas_adsorb_results', **my_gcmc_props)
```

The method `css.run()` will trigger the simulations from run queue. In our case, it will be one GCMC simulation.


### Working with simulation results

The results of Cassandra simulations are written to the **cassandra.system** field that is the object of **pysimm.system** type. 
The **pysimm.system** has few methods to write itself to the text-formatted files (such as .lmps or .xyz) for further 
molecular simulation/visualization treating.

```python
css.system.write_lammps('gas_adsorb.lmps')
css.system.write_xyz('gas_adsorb.xyz')
```
