Examples
========

These examples demonstrate some of the most useful features of pysimm, including structure generation, force field parameter assignment, and LAMMPS simulation. Each example consists of a python script and a README markdown file that explains each line of code. The first 6 examples can be performed using either the Dreiding or GAFF2 force fields. The GAFF2 examples require an extra, optional LAMMPS package, "user-misc", and will not work with versions of LAMMPS without this package installed. Users are encouraged to start with the Dreiding examples, and if desired, install a working version of LAMMPS with the "user-misc" package installed using the complete_install script supplied with the pysimm source code.

Dreiding Example Requirements:
* LAMMPS (molecule, kspace packages)
* Dreiding force field data file (pysimm/dat/forcefields/dreiding.xml)
* Internet access

GAFF2 Example Requirements:
* LAMMPS (molecule, kspace packages, user-misc)
* GAFF2 force field data file (pysimm/dat/forcefields/gaff2.json)
* Internet access