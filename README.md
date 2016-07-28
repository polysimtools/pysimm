pysimm
======

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

pysimm stands for <b>Py</b>thon <b>S</b>imulation <b>I</b>nterface for <b>M</b>olecular <b>M</b>odeling. It is an open-source python package designed to assist in the setup and executation of molecular simulations through high level APIs and abstraction from underlying third party simulation software.

Getting Started
===============

To get started, clone the repository, cd into the new directory, and install using python setuptools:

```
git clone https://github.com/polysimtools/pysimm
cd pysimm
sudo python setup.py install
```

This adds the pysimm package to your PYTHONPATH, and adds a command line tool to your PATH. 

Integration with LAMMPS
=======================

pysimm can integrate seamlessly with parts of the LAMMPS simulation software package through the pysimm.lmps module. To configure the integration, locate your LAMMPS executable. If the path to your LAMMPS executable is /usr/bin/lmp_mpi, add this path as an environment variable "LAMMPS_EXEC":

```echo "export LAMMPS_EXEC=/usr/bin/lmp_mpi" >> ~/.bashrc```

If the LAMMPS_EXEC environment variable is not set, you will see a warning when importing the lmps module. Once configured, try the examples in the repository which highlight some of the features of pysimm.

Complete Installation (pysimm and LAMMPS)
=========================================

<i>NOTE: The complete_install.py script is designed to work on debian-based linux machines only. It uses apt-get to install dependencies for pysimm and LAMMPS.</i>

Included in the repository is a python script complete_install.py that will configure pysimm, install LAMMPS from their git repository, and configure the integration between the two pieces of software. First clone the pysimm repository, and run complete_install.py like such:

```python complete_install.py --pysimm /path/to/pysimm/repo --lammps /path/to/new/lammps/directory```

Afterwords be sure to source your ~/.bashrc file:

```source ~/.bashrc```