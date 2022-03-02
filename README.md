pysimm
======

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

pysimm stands for <b>Py</b>thon <b>S</b>imulation <b>I</b>nterface for <b>M</b>olecular <b>M</b>odeling. It is an open-source python package designed to assist in the setup and executation of molecular simulations through high level APIs and abstraction from underlying third party simulation software.

Documentation for the python package can be found at http://pysimm.org/documentation

### Sections
1. [Getting Started](#getting-started)
2. [Integration with LAMMPS](#integration-with-lammps) 
3. [Integration with Cassandra](#integration-with-cassandra)
4. [Complete Installation](#complete-installation-pysimm-and-lammps)
5. [Using Docker image](#using-docker-image)
6. [Acknowledgments](#acknowledgments)


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

Correct configuration can be checked by using the following code:

```
from pysimm import lmps
lmps.check_lmps_exec()
```

If the LAMMPS_EXEC environment variable is not set, you will see a warning. Once configured, try the examples in the repository which highlight some of the features of pysimm.


Integration with Cassandra 
==========================
<a name="CSintegration"></a>
Cassandra is a versatile software for the Monte-Carlo simulations suitable for the variety of the molecular systems. For more details please refer the [web-page](https://cassandra.nd.edu/) of the Cassandra project. This explains how to integrate already installed Cassandra software to the pysimm. The integration is done through the environment variable CASSANDRA_EXEC. It defines the path to the Cassandra executable (it might be either serial or the parallel Cassandra compilation). If, for example, the Cassandra executable named *cs_omp_gfort.exe* was installed to */usr/lib/cassandra/*, the following should be set:

```export CASSANDRA_EXEC=/usr/lib/cassandra/Src/cs_omp_gfort.exe```

The cassandra module of the pysimm will check whether the CASSANDRA_EXEC is set right before the running a simulation. Additionally, in the same way as in the lammps module, the setting can be checked using the code:

```
from pysimm import cassandra
cassandra.check_cs_exec()
```
<i>NOTE: The parallel Cassandra version uses OpenMP, so the number of parallel tasks (number of threads) is normally controlled through OMP_NUM_THREADS environment variable. Set it to required value if you want to use the specific number of parallel threads.</i>


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


Using Docker image
==================

The root folder now contains a 'Dockerfile' that will help you to create a Docker 
image of Debian 10 with pysimm and pre-installed LAMMPS, 
CASSANDRA, Zeo++ v0.3, and PoreBlazer v4.0.  

To compile the Docker image from the file run the following from the root pysimm directory 
(optionally changing 'my_tag' to any text tag you like):

```commandline
 docker build -t pysimm:my_tag -f Dockerfile .
```

Please refer the [Docker documentation](https://docs.docker.com/engine/reference/commandline/build/) 
for the detailed description of the `build` function.  

If the build is successful the list of your Docker images will contain freshly built 
pysimm image. The full list can be shown by:
```commandline
 docker images
```

Finally, to run the corresponding pysimm image call:  
```commandline
 docker run -it pysimm:my_tag bash
```
Please see [Docker reference page](https://docs.docker.com/engine/reference/run/) for the 
detailed description of `docker run` command.

The pysimm source files are kept in `/usr/local/pysimm` folder. Thus you can quickly test the 
LAMMPS or CASSANDRA modules by running one of the examples.
  
```commandline
 cd /usr/local/pysimm/Examples/08_ethanol_acetone_mixture
 python run.py
```

To transfer output files from the docker image back to host one can use `docker cp` command. 
```commandline
 docker cp ContainerName:/container/path/to/file /host/path/to/file
```

and the `ContainerName` one can use either the name or ID of a running container 
which can be listed by calling `docker ps`.
 

Acknowledgments
================
This material is based upon work supported by: 
 * National Science Foundation under Grant No. (ACI-1613155)
 * U.S. Departement of Energy under Grant No. (DE-FG02-17ER16362) 
