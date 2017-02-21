pysimm
======

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

pysimm stands for <b>Py</b>thon <b>S</b>imulation <b>I</b>nterface for <b>M</b>olecular <b>M</b>odeling. It is an open-source python package designed to assist in the setup and executation of molecular simulations through high level APIs and abstraction from underlying third party simulation software.

Documentation for the python package can be found at http://pysimm.org/documentation

Getting Started
===============

Those who need to install LAMMPS (or wish to build a new version of LAMMPS with the optional packages utilized in the examples included in the repository) should skip to the Complete Installation section below, which provides instructions on how to setup pysimm, build a LAMMPS version suitable for the examples included here, and configure the integration between the two pieces of software. If you have already installed LAMMPS and simply wish to setup pysimm, continue following the directions below.

To get started, clone the repository, cd into the new directory, and install using complete_install.py. The --pysimm command line argument passed to the script should be the directory in which you cloned the pysimm repository (one directory up). The following example assumes you cloned the repository in your home directory.

```
git clone https://github.com/polysimtools/pysimm
python pysimm/complete_install.py --pysimm $PWD
```

This adds the pysimm package to your PYTHONPATH, and adds a pysimm command line tool to your PATH. Parts of pysimm require the use of the numpy package. To use the complete_install script to install numpy as well, include the --apt-install command line argument.

Integration with LAMMPS
=======================

If you are using your own build of LAMMPS, be sure that the following packages were included in your installation as some functionality in the example scripts require some subset of these packages:
  -  molecule
  -  class2
  -  kspace
  -  user-misc
  -  misc
  -  qeq
  -  manybody

pysimm can integrate seamlessly with parts of the LAMMPS simulation software package through the pysimm.lmps module. To configure the integration, locate your LAMMPS executable. For example, if the path to your LAMMPS executable is /usr/bin/lmp_mpi, add this path as an environment variable "LAMMPS_EXEC":

```echo "export LAMMPS_EXEC=/usr/bin/lmp_mpi" >> ~/.bashrc```

If the LAMMPS_EXEC environment variable is not set, you will see a warning when importing the lmps module. Once configured, try the examples in the repository which highlight some of the features of pysimm.

Complete Installation (pysimm and LAMMPS)
=========================================

*** Please be aware LAMMPS is a separate software, developed separately, and protected under a separate license. The following information is installation directions <i>suggested</i> by pysimm developers that may be used to generate a functional version of LAMMPS for integration with pysimm. LAMMPS has much more functionality than is represented by these installation directions, and those interested in these functionalities should visit [LAMMPS documentation](http://lammps.sandia.gov/doc/Manual.html) to learn more. 

<i>NOTE: The complete_install.py script is designed to work on debian-based linux machines only. It uses apt-get to install dependencies for pysimm and LAMMPS.</i>

Included in the repository is a python script complete_install.py that will configure pysimm, install LAMMPS from their git repository, and configure the integration between the two pieces of software. First clone the pysimm repository, and run complete_install.py. You must provide the path prefix to the recently cloned pysimm directory and the prefix for the new lammps source code directory. The following assumes pysimm was cloned into your home directory and will also install lammps in your home directory:

```
git clone https://github.com/polysimtools/pysimm
python pysimm/complete_install.py --pysimm $PWD --lammps $PWD
```

Afterwords be sure to source your ~/.bashrc file:

```source ~/.bashrc```