import os
import sys
import errno
import argparse
from subprocess import call, PIPE, Popen

HOME_DIR = os.environ.get('HOME')


def install_pysimm(prefix):
    os.chdir(prefix)
    if os.path.isdir('pysimm'):
        print('pysimm directory already exists...assuming it is the pysimm repository and continuing...')
    else:
        call('git clone https://github.com/polysimtools/pysimm', shell=True)
    if not os.path.isfile('pysimm/complete_install.py'):
        print('assumption about pysimm repository existing was wrong; exiting...')
        exit()
    call("echo export PYTHONPATH='$PYTHONPATH':{} >> {}".format(os.path.join(prefix, 'pysimm'),
                                                                os.path.join(HOME_DIR, '.bashrc')),
         shell=True)
    call("echo export PATH='$PATH':{} >> {}".format(os.path.join(prefix, 'pysimm', 'bin'),
                                                    os.path.join(HOME_DIR, '.bashrc')),
         shell=True)


def apt_update():
    call('apt-get update', shell=True)


def apt_install(*packages):
    call('apt-get -y install {}'.format(' '.join(packages)),
         shell=True)


def install_lammps(prefix, *packages):
    os.chdir(prefix)
    call('git clone -b unstable https://github.com/lammps/lammps.git lammps', shell=True)
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
         
def install_ambertools(dir_):
    os.chdir(dir_)
    call("echo export AMBERHOME={} >> {}".format(dir_, os.path.join(HOME_DIR,'.bashrc')),
         shell=True)
    os.environ['AMBERHOME'] = dir_
    call('./configure gnu', shell=True)
    call('make install', shell=True)
    call("echo export ANTECHAMBER_EXEC={} >> {}".format(os.path.join(dir_, 'bin', 'antechamber'),
                                                   os.path.join(HOME_DIR,'.bashrc')),
         shell=True)
         
def install_openbabel():
    apt_install('libopenbabel4', 'libopenbabel-dev', 'openbabel', 'python-openbabel')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--apt-update', dest='apt_update', action='store_true', default=False)
    parser.add_argument('--apt-install', dest='apt_install', action='store_true', default=False)
    parser.add_argument('--pysimm', dest='pysimm_prefix', default=HOME_DIR)
    parser.add_argument('--lammps', dest='lammps_prefix', default=None)
    parser.add_argument('--lammps-packages', dest='lammps_packages', nargs='*',
                        default=['molecule', 'class2', 'kspace', 'user-misc', 'misc', 'qeq', 'manybody'])
    parser.add_argument('--amber-tools', dest='ambertools_dir', default=None)
    parser.add_argument('--openbabel', dest='openbabel', action='store_true', default=False)
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

    if bool(args.apt_update):
        apt_update()

    if args.pysimm_prefix:
        if args.apt_install:
            apt_install('git', 'python-numpy', 'python-matplotlib')
        mkdir_p(args.pysimm_prefix)
        install_pysimm(args.pysimm_prefix)

    if args.lammps_prefix:
        if args.apt_install:
            apt_install('make git g++', 'libopenmpi-dev', 'openmpi-bin')
        mkdir_p(args.lammps_prefix)
        install_lammps(args.lammps_prefix, *args.lammps_packages)
        
    if args.ambertools_dir:
        if args.apt_install:
            apt_install('make', 'csh', 'gfortran', 'libopenmpi-dev', 'openmpi-bin', 'xorg-dev', 'xserver-xorg')
        install_ambertools(args.ambertools_dir)
        
    if args.openbabel:
        install_openbabel()
        
    os.chdir(HOME_DIR)

