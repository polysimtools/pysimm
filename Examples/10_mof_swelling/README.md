Example 10: Hybrid Monte Carlo – Molecular Dynamics (MC-MD) simulations with PySIMM
========================================================================================
by Alexander Demidov and Michael E. Fortunato

### General description of the hybrid MC-MD simulation. Modules import

The example describes the work with the pysimm app that implements the hybrid MC-MD approach. Under the "hybrid MC-MD approach" here we understand the following serial iterative set of simulations. First, the Monte-Carlo (MC) simulation is performed, inserting the gas molecules into the simulation system with constant volume temperature and chemical potential of inserted species (&#956;VT). This is followed by the molecular dynamics (MD) simulation step at constant pressure, temperature and number of particles (NPT). After the MD is finished, the next step of the &#956;VT MC with the new size of the simulation box is performed, and so on.

In the pySIMM the hybrid MC-MD simulations are implemented in **apps.mc_md** application. The application in turn communicates with [Cassandra](https://cassandra.nd.edu) and [LAMMPS](http://lammps.sandia.gov)  programs through **pysimm.cassandra** and **pysimm.lmps** modules correspondingly automatically setting up the required for simulation data.

```python
from pysimm.apps import mc_md
from pysimm import system
import os
```

The example uses the hybrid MC-MD for simulation of swelling of an IRMOF-1 metal-organic framework with pure methane. 

### Simulation system setup

To start simulations, the application should get two **pysimm.system** objects. First, *frame*, is any molecular structure that presumed to be fixed during the MC simulations, and the second, *gas1*, that represents single gas molecule to be inserted by the MC. 

```python
frame = system.read_lammps(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'irmof1_drei.lmps'))
frame.forcefield = 'dreiding-lj'
gas1 = system.read_lammps(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ch4.lmps'))
gas1.forcefield = 'trappe/amber'
```

### Simulation properties setup

Additionally the application requires two dictionaries *mc_props* and *md_props* that describe the Cassandra and LAMMPS simulation settings correspondingly. The keywords of the property dictionaries correspond to the keywords provided to the Cassandra and LAMMPS programs directely.

```python
mc_props = {'rigid_type': False,
            'max_ins': 1000,
            'Chemical_Potential_Info': -30.09,
            'Temperature_Info': 300,
            'Run_Type': {'steps': 250},
            'CBMC_Info': {'rcut_cbmc': 2.0},
            'Simulation_Length_Info': {'run': 10000,
                                       'coord_freq': 10000,
                                       'prop_freq': 500},
            'VDW_Style': {'cut_val': 9.0},
            'Charge_Style': {'cut_val': 9.0},
            'Property_Info': {'prop1': 'energy_total',
                              'prop2': 'pressure'}}

md_props = {'temp': 300,
            'pressure': {'start': 15,
                         'iso': 'iso'},
            'timestep': 0.25,
            'cutoff': 9.0,
            'length': 10000,
            'thermo': 2500,
            'dump': 2500,
            'print_to_screen': False}
```

### Running the application

The application itself has two settings that are provided in the form of keyword arguments: the number of iterative loops *mcmd_niter* (default is 10) and the relative path to all simulation results *sim_folder* (default is *results*).  The call of **mc_md()** function will automatically start the MC-MD simulations.

```python
sim_result = mc_md.mc_md(gas1, polymer, mcmd_niter=3, sim_folder='results',  mc_props=mc_props, md_props=md_props)
```

The application returns the **pysimm.system** object that is a dump of simulated system after the final iteration.