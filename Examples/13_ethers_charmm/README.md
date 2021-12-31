Example 13: Water solvation of a linear ether using CHARMM forcefield
========================================================================================================================  
by Alexander Demidov
    
### Summary  

The example shows a simple workflow that constructs a simulation system with a single solute molecule (ethyl-propyl-ether) 
submerged into an explicit solvent (water). The ether molecule is typed with CHARMM forcefield and the PySIMM system is 
prepared such that it is ready for simulation. Short MD simulations are set up in this example with LAMMPS.

      
### 1. Read data
Initially, ethyl-propyl-ether molecule read from a .mol file which contains information only about atoms (as elements) and 
bonds between atoms. Let us assign al bonds order of 1 (it is a simple ether molecule)

```python
for b in sst.bonds:
    b.order = 1

```   
Element names of all atoms bonds and bond orders is the information needed for the force field automatic typing in PySIMM. 
 In the example, we use CHARMM force field and simple Gasteiger algorithm to automatically assign atom types
 
 ```python
ff = forcefield.Charmm()
sst.apply_forcefield(ff, charges='gasteiger')
sst.set_charge()
```

Additionally, we expand the simulation box to 2.5 nm to make its linear size longer than twice the cutoff distance of the force 
field (1.2 nm). Next, read in the solvent molecule (water, cTIP3P model).  


### 2. Combine simulation systems

Simplest way to construct the initial position of water molecules would be to put them into nodes of the regular grid. 
To calculate how many cells the grid should contain let's assume that density of initial system should be 1 
(gram / cm<sup>3</sup>).
```python
ngrid = numpy.floor(((bxSize ** 3) * 0.6022 / solvnt.mass) ** (1.0 / 3.0))
```

Finally, let's put all water molecules in their initial positions copying them from previously read system. 
Note, that we will not add a water molecule if its center closer than 1.7 A to any of the ether atoms. 

```python
rng = numpy.linspace(sst.dim.xlo, sst.dim.xhi, int(ngrid) + 1)
count = 0
for p in rng[:-1]:
    for q in rng[:-1]:
        for t in rng[:-1]:
            flags = []
            for prt in sst.particles:
                dist = numpy.linalg.norm(numpy.array([prt.x, prt.y, prt.z]) - numpy.array([p, q, t]))
                flags.append(dist > 1.7)

            if all(flags):
                count += 1
                tmp = solvnt.copy(dx=p, dy=q, dz=t)
                sst.add(tmp, change_dim=False, update_properties=False)
```

Because CHARMM normally works with explicit Lennard-Jones parameters for different types of atoms rather than utilizing 
mixing rules, one need to manually update non-diagonal LJ parameters of the system:

```python
ff.assign_extra_ljtypes(sst)
``` 

### 3. Setup LAMMPS simulations
In general, the setup of LAMMPS simulation is as usual and includes initialization of simulation, adjustment of output 
settings, and setup of molecular dynamic settings. However, two additional extra modifications are important for this example:   

**i)** In LAMMPS aforementioned explicit Lennard-Jones parameters for different types of atoms are defined with the `pair_coeff` 
command that is written explicitly into the input file. 
This can be concisely done with pysimm for the system of any number of those additional parameters:

```python
for nd_lj in sst.nondiag_lj_types:
    sim.add_custom('pair_coeff {} {} {}'.format(' '.join(map(str, nd_lj.atm_types)), nd_lj.epsilon, nd_lj.sigma))

```  
**ii)** Because we work with the rigid model of water we want to set up and additional shake fix, which will fix the 
length of H-O bond and H-O-H angle. 

```python
sim.add_custom('fix shck_fix all shake 0.0001 10 0 b {:} a {:}\n'.format(sst.bond_types.get('H,O')[0].tag,
                                                                        sst.angle_types.get('H,O,H')[0].tag))
``` 

Finally, we run the simulations and save the snapshot of the equilibrated system in the format of the .lmps file .
