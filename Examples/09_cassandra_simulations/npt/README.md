Example 9.2: The Isothermal-Isobaric Monte Carlo simulations using Cassandra
============================================================================
by Alexander Demidov and Michael E. Fortunato

### Importing pysimm modules/packages

The example illustrates the basic principles of utilizing the **pysimm.cassandra** module. The code sets up the Isothermal-Isobaric (NPT) Monte Carlo simulations to be performed by the [Cassandra](https://cassandra.nd.edu) software. The example shows how to change Cassandra simulation settings and save the simulation results.

First, let's import the required pysimm modules: system and cassandra

```python
from pysimm import system, cassandra
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Creating the system for simulations

The cassandra module of the pysimm works with **pysimm.system** objet modyfing it by the results of Cassandra simulations. For this example the initial system can be an empty cubic cell with the size of 60 &#8491; cenetred at the origin of the coordinate system. Additionally we explicitly specify the forcefield model for all particles that will be inserted into the system.

```python
sst = system.System()
sst.dim = system.Dimension(dx=45, dy=45, dz=45, center=[0, 0, 0])
sst.forcefield = 'trappe/amber'
```

The setup of a Monte Carlo simulations should include also the information about entities to be modified by the Monte Carlo algorithm. This example works with  ethyne (CH<sub>2</sub>--CH<sub>2</sub>) molecule using TraPPE-UA forcefield. For the **pysimm.cassandra** module they should be provided in the form of the **pysimm.system** objects.

The example folder contains the LAMMPS files that describe the ga molecule, so the corresponding  **pysimm.system** objects can be created by reading the files by the **system.read_lammps()** method.

```python
molec = system.read_lammps('c2h4.lmps')
molec.forcefield = 'trappe/amber'
```

The system object is passed to the constructor of the main cassandra simulation object 

```python
css = cassandra.Cassandra(sst)
```

### Creating/changing the simulation properties

The properties of the Cassandra simulation for the pysimm.cassandra module are provided in the form of a dictionary in which the key is the property name and the value (depending on the property itself) can be integer, float, string or another dictionary. Thus, all required Cassandra simulation properties can be specified directly in the Python script file by creating the dictionary. Alternatively, pysimm.cassandra object can read and interpret the proper Cassandra input file. If the file does not specify all properties required for certain type of simulations, they with set to the default ones.

Let's read some data from the **props.inp** file that is in the same example directory

```python
my_gcmc_props = css.read_input('props.inp')
```

Let's explicitely change some properties that will define the NPT simulations directly: pressure, initial configuration, simulation length, and output style:
```python
npt_props['Pressure_Info'] = 25  # Simulated pressure in bars
npt_props['Start_Type'] = {'start_type': 'make_config', 'species': 300}
npt_props['Simulation_Length_Info'] = {'run': 10000}
npt_props['Property_Info'] = {'prop1': 'energy_total',
                              'prop2': 'volume',
                              'prop3': 'mass_density'}
```

### Setup and run the simulations

The set-up of NPT simulations in **pysimm.cassandra** module are implemented in NPT class. The way to add it to the simulation queue is to call **css.add_npt_mc()** method. The method has few keyword arguments:
* *species*: the molecules to be inserted;
* *out_folder*: relative path for all Cassandra files such as log files, .xyz and .chk files

```python
cs.add_npt_mc(species=molec, is_rigid=True, out_folder='results', **npt_props)
```

The method `css.run()` will trigger the simulations from run queue. In our case it will be one GCMC simulation.


### Working with simulation results

The results of Cassandra simulations are written to the **cassandra.system** field that is object of **pysimm.system** type.  The **pysimm.system** has few methods to write itself to the text-formatted file (e.g. such as .lmps or .xyz).


```python
cs.system.write_lammps('final_conf.lmps')
```
