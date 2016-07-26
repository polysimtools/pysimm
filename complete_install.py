import os
import sys
import errno
import argparse
from subprocess import call, PIPE, Popen

HOME_DIR = os.environ.get('HOME')


def install_pysimm(prefix):
    os.chdir(prefix)
    call('git clone https://github.com/polysimtools/pysimm', shell=True)
    call("echo export PYTHONPATH='$PYTHONPATH':{} >> {}".format(prefix,
                                                                os.path.join(HOME_DIR, '.bashrc')),
         shell=True)
    call("echo export PATH='$PATH':{} >> {}".format(os.path.join(prefix, 'bin'),
                                                    os.path.join(HOME_DIR, '.bashrc')),
         shell=True)


def apt_update():
    call('sudo apt-get update', shell=True)


def apt_install(*packages):
    call('sudo apt-get -y install {}'.format(' '.join(packages)),
         shell=True)


def install_lammps(prefix, *packages):
    os.chdir(prefix)
    call('git clone git://git.lammps.org/lammps-ro.git lammps'.format(prefix), shell=True)
    os.chdir(os.path.join(prefix,'lammps','src'))
    for package in packages:
        call('make yes-{}'.format(package), shell=True)
    call('make mpi', shell=True)
    call("echo export PATH='$PATH':{} >> {}".format(os.path.join(prefix, 'lammps', 'src'),
                                                    os.path.join(HOME_DIR,'.bashrc')),
         shell=True)
    call("echo export LAMMPS_EXEC={} >> {}".format(os.path.join(prefix, 'lammps', 'src', 'lmp_mpi'),
                                                   os.path.join(HOME_DIR,'.bashrc')),
         shell=True)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pysimm', dest='pysimm_prefix', default=HOME_DIR)
    parser.add_argument('--lammps', dest='lammps_prefix', default=None)
    parser.add_argument('--lammps-packages', dest='lammps_packages', nargs='*',
                        default=['molecule', 'class2', 'kspace'])
    return parser.parse_args()


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


if __name__ == '__main__':

    args = parse_args()

    apt_update()

    if args.pysimm_prefix:
        apt_install('git', 'python-numpy', 'python-matplotlib')
        mkdir_p(args.pysimm_prefix)
        install_pysimm(args.pysimm_prefix)

    if args.lammps_prefix:
        apt_install('make git g++', 'libopenmpi-dev', 'openmpi-bin')
        mkdir_p(args.lammps_prefix)
        install_lammps(args.lammps_prefix, *args.lammps_packages)
    os.chdir(HOME_DIR)
