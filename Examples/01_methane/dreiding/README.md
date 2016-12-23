Example 1: Creating Methane Atom-by-Atom
========================================
by Michael E. Fortunato

### Importing pysimm modules/packages

This example code creates a new molecular system from scratch and builds a methane molecule atom-by-atom. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield.

`from pysimm import system, lmps, forcefield`

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Preparing an empty system

First we will create an empty **system.System** object and store this in a variable **s**. This system object will contain and organize all of our molecular data:

`s = system.System()`

Next we need to add a molecule to our system. By default, our system **s** has a container object **s.molecules** that we need to add our new molecule to. We create a new molecule object, **system.Molecule()**, and pass this to the molecule container class method **s.molecules.add()**. This function returns the newly added object to the container, and we store this in variable **m**.

`m = s.molecules.add(system.Molecule())`

### Obtaining force field parameters

Now that we have a place to contain all of our molecular data, we need to obtain data describing interactions between atoms from a forcefield. A **forcefield.Forcefield** object contains the parameters necessary to calculate the energy of your system as well as logical typing rules to assign atom types to particles in our system. In this example, we will use the Dreiding force field and instantiate a Dreiding force field object and store it in variable **f**:

`f = forcefield.Dreiding()`

This example assumes the user already knows which Dreiding atom types are required for a methane molecule: "C_3" and "H_". Let's get **system.ParticleType** objects from our **forcefield.Dreiding** object that represent these atom types. Our **forcefield.Dreiding** object **f** has a container object **f.particle_types** that organizes **system.ParticleType** objects and provides a method **f.particle_types.get()** to retrieve specific types based on names. The function returns a list of **system.ParticleType** objects, but in this example we know only one object is returned so we access the first element of the list. Instead of adding the **system.ParticleType** from **f**, we create a copy, and pass this newly created **system.ParticleType** object to a method that adds it to the particle type container object in our system **s.particle_types**. The **s.particle_types.add()** method returns the newly added object, which we store as the Dreiding atom type object representation.

```
dreiding_C_ = s.particle_types.add(f.particle_types.get('C_3')[0].copy())
dreiding_H_ = s.particle_types.add(f.particle_types.get('H_')[0].copy())
```

### Adding atoms to a system/molecule

Now we have a System object, **s**, a Molecule object that is stored in our system, **m**, and two Dreiding atom type objects, **dreiding_C_** and **dreiding_H_**. Let's start adding atoms.

First we create the carbon atom (particle) at the origin. For now we add the particle with zero charge, and we will derive Gasteiger charges later. We instantiate a **system.Particle** object using keyword arguments to set x, y, and z coordinates, the charge, the molecule this particle is part of, and most importantly, the particle type, **dreiding_C_**. Notice this is a reference to our **system.ParticleType** object that is in our **system.System** object, not the **forcefield.Dreiding** object. Again we use a container add method to add our new particle to the particles container, and store the newly added particle in variable **c1**.

`c1 = s.particles.add(system.Particle(type=dreiding_C_, x=0, y=0, z=0, charge=0, molecule=m))`

The **system.System** class has a method that allows a user to add a new particle that should be bonded to an existing particle in the system **system.System.add_particle_bonded_to**. If a force field object is also passed to this method, new bonds, angles, and dihedrals will be created as necessary. The location of the new atom is random within a carbon bond length radius around the existing particle, and we will later use LAMMPS to perform structural optimization. For each hydrogen atom we want to add, we create a new **system.Particle** object, with zero charge and the **dreiding_H_** type. We also need to identify this will be part of the molecule that already exists, **m**. This **system.Particle** object is passed to **s.add_particle_bonded_to()** as the first argument. The second argument is the particle that already exists in our system, in this case **c1**. The third argument is a reference to our **forcefield.Dreiding** object. The **s.add_particle_bonded_to()** method returns the newly added object, which we store as **h1**, **h2**, **h3** and **h4**.

```
h1 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
h2 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
h3 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
h4 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
```

Lastly we can derive Gasteiger charges using the **system.System.apply_charges()** class method, providing a **forcefield.Forcefield** object and identifying we would like to use the Gasteiger approach.

`s.apply_charges(f, charges='gasteiger')`

### Setting up a simulation box and optimizing the structure

We've added particles to our system but there is currently no simulation box. The **system.System** class has a method to construct a simulation box surrounding the atoms it contains, and optionally we can add a padding around our molecule. Here we opt to use a padding of 10 Angstroms.

`s.set_box(padding=10)`

Before we optimize our structure, LAMMPS will need to know what type of pair, bond, and angle interactions we are using. These are defined as attributes of our **system.System** object.

```
s.pair_style='lj'
s.bond_style='harmonic'
s.angle_style='harmonic'
```

We'll use the fire minimization algorithm implemented in LAMMPS. The **lmps** module in pysimm contains convenience methods to configure and execute a simulation. In this case we will use **lmps.quick_min()** passing our system object **s**, the min_style we want to use, and give each simulation a name so that the log files have identifiable names.

```
lmps.quick_min(s, min_style='fire', name='fire_min')
```

### Writing out a system to various file formats

The **system.System** class has various method to format our system data into different file types. Here we write xyz, yaml, lammps data, and chemdoodle json file formats.

```
s.write_xyz('methane.xyz')
s.write_yaml('methane.yaml')
s.write_lammps('methane.lmps')
s.write_chemdoodle_json('methane.json')
```