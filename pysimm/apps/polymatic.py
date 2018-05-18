# polymatic module; calls Polymatic perl code and LAMMPS

import os
import sys
import shlex
import shutil
from time import strftime
from subprocess import Popen, PIPE

from pysimm import system, lmps

rappture = True
try:
    import Rappture
except ImportError:
    rappture = False


def pack(script, file_in, nrep, boxl, file_out):
    """pysimm.apps.polymatic.pack

    Calls Polymatic random packing code

    Args:
        script: name of packing script
        file_in: list of file names of reference molecules to pack
        nrep: list of number of monomers for each reference molecule
        boxl: length of one dimension of simulation box for random packing
        file_out: name of output file (packed system)

    Returns:
        output from perl code
    """
    if not isinstance(file_in, list):
        file_in = [file_in]
    if not isinstance(nrep, list):
        nrep = [nrep]
    if len(file_in) != len(nrep) or len(file_in) == 0:
        return False

    cmd = 'perl %s -i ' % script
    cmd += '%s ' % str(len(file_in))
    for i in range(len(file_in)):
        cmd += '%s %s ' % (file_in[i], nrep[i])
    cmd += '-l %s -o %s' % (boxl, file_out)

    o, e = Popen(shlex.split(cmd),
                 stdin=PIPE,
                 stdout=PIPE,
                 stderr=PIPE).communicate()

    if not e and o:
        return o
    else:
        return False


def polymatic(script, file_in, file_out):
    """pysimm.apps.polymatic.polymatic

    Calls Polymatic code. polym.in and types.txt are assumed to exist.

    Args:
        script: name of Polymatic script
        file_in: initial system file name
        file_out: final system file name

    Returns:
        output from perl code
    """
    cmd = ('perl %s -i %s -s polym.in -t types.txt -o %s'
           % (script, file_in, file_out))
    o, e = Popen(shlex.split(cmd),
                 stdin=PIPE,
                 stdout=PIPE,
                 stderr=PIPE).communicate()

    if not e and o and o.split()[0] is not 'Error:':
        return True
    elif not e and o:
        return o
    else:
        return False


def run(settings):
    """pysimm.apps.polymatic.run

    Runs Polymatic algorithm.

    Args:
        settings: object containing Polymatic settings

    Returns:
        (True/False, :class:`~pysimm.system.System`)
    """
    if rappture:
        Rappture.Utils.progress(0, 'Initializing Polymatic...')
    bonds = 0

    os.mkdir('logs')

    polymatic(os.path.join(settings.polymatic_dir, 'polym_init.pl'),
              'data.lmps',
              'step_000.lmps')

    s = system.read_lammps('step_000.lmps', quiet=True)
    s.read_type_names('types.txt')
    s.write_lammps('temp.lmps')

    if rappture:
        Rappture.Utils.progress(0, '%s/%s bonds made: '
                                   'optimizing initial structure...'
                                % (bonds, settings.polym.max_bonds))

    if not lmps_min(s, 'initial optimization', settings):
        s.write_lammps('temp.lmps')
        polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
                  'temp.lmps',
                  'final.lmps')
        return False, s

    s.write_lammps('step_000.lmps')
    s.write_lammps('temp.lmps')

    while bonds < settings.polym.max_bonds:
        attempt = 0
        while not polymatic(os.path.join(settings.polymatic_dir, 'polym.pl'),
                            'temp.lmps',
                            'temp.lmps'):
            attempt += 1
            if rappture:
                Rappture.Utils.progress(int(float(bonds)/settings.
                                            polym.max_bonds *
                                            100),
                                        '%s/%s bonds made: attempt #%s to make '
                                        'new bond'
                                        % (bonds, settings.polym.max_bonds,
                                           attempt))
            s = system.read_lammps('temp.lmps', quiet=True)
            s.read_type_names('types.txt')

            if not lmps_step_md(s, bonds, attempt, settings):
                s.write_lammps('temp.lmps')
                polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
                          'temp.lmps',
                          'final.lmps')
                return False, s
            s.write_lammps('temp.lmps')

            if attempt >= settings.polym.max_md:
                break

        if attempt >= settings.polym.max_md:
            break

        bonds += 1

        if rappture:
            Rappture.Utils.progress(int(float(bonds)/settings.polym.max_bonds
                                        * 100),
                                    '%s/%s bonds made: '
                                    'optimizing newly formed bond'
                                    % (bonds, settings.polym.max_bonds))

        s = system.read_lammps('temp.lmps', quiet=True)
        s.read_type_names('types.txt')

        print('%s: bond %s made successfully' % (strftime('%H:%M:%S'), bonds))
        sys.stdout.flush()
        if not lmps_min(s, 'bond %s optimization' % bonds, settings):
            s.write_lammps('temp.lmps')
            polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
                      'temp.lmps',
                      'final.lmps')
            return False, s
        s.write_lammps('step_%03d.lmps' % bonds)
        s.write_lammps('temp.lmps')

        if (bonds % settings.polym.cycle == 0 and
                (bonds / settings.polym.cycle) % settings.polym.npt_freq == 0):
            if rappture:
                Rappture.Utils.progress(int(float(bonds)/settings.
                                            polym.max_bonds
                                            * 100),
                                        '%s/%s bonds made: '
                                        'performing npt cycle md'
                                        % (bonds, settings.polym.max_bonds))
            if not lmps_cycle_npt_md(s, bonds, settings):
                s.write_lammps('temp.lmps')
                polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
                          'temp.lmps',
                          'final.lmps')
                return False, s
            s.write_lammps('temp.lmps')

        elif bonds % settings.polym.cycle == 0:
            if rappture:
                Rappture.Utils.progress(int(float(bonds)/settings.
                                            polym.max_bonds
                                            * 100),
                                        '%s/%s bonds made: '
                                        'performing nvt cycle md'
                                        % (bonds, settings.polym.max_bonds))
            if not lmps_cycle_nvt_md(s, bonds, settings):
                s.write_lammps('temp.lmps')
                polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
                          'temp.lmps',
                          'final.lmps')
                return False, s
            s.write_lammps('temp.lmps')

    if rappture:
        Rappture.Utils.progress(99, 'Finalizing Polymatic')

    polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
              'temp.lmps',
              'final.lmps')

    return True, s


def lmps_min(s, name, settings):
    """pysimm.apps.polymatic.lmps_min

    Runs LAMMPS minimization for the Polymatic algorithm.

    Args:
        s: :class:`~pysimm.system.System` to minimize
        name: name of simulation
        settings: object containing Polymatic settings

    Returns:
        result from :func:`~pysimm.lmps.minimize`
    """
    if settings.polym.min.cluster:
        nanohub = {'cores': int(settings.polym.min.nanohub_cores),
                   'walltime': int(settings.polym.min.nanohub_walltime)}
        log_name = '%s' % '_'.join(name.split())
    else:
        nanohub = {}
        log_name = 'logs/%s' % '_'.join(name.split())

    if settings.polym.min.user_input:
        sim = lmps.Simulation(s, name='initial optimization',
                 print_to_screen=False, nanohub=nanohub, custom=True)
        sim.add(settings.polym.min.min_in)
        sim.run(np=settings.np, nanohub=nanohub)
    else:
        sim = lmps.Simulation(s, name='initial optimization',
            print_to_screen=False, nanohub=nanohub,
            log=log_name
        )
        sim.add(lmps.Init(cutoff=settings.polym.min.nb_cutoff, forcefield=settings.forcefield))
        sim.add_min(
            min_style='sd',
            etol=settings.polym.min.sd_etol,
            ftol=settings.polym.min.sd_ftol,
            maxiter=settings.polym.min.sd_maxiter,
            maxeval=settings.polym.min.sd_maxeval,
        )
        sim.add_min(
            min_style='cg',
            etol=settings.polym.min.cg_etol,
            ftol=settings.polym.min.cg_ftol,
            maxiter=settings.polym.min.cg_maxiter,
            maxeval=settings.polym.min.cg_maxeval,
        )
        sim.run(np=settings.np, nanohub=nanohub)

    if settings.polym.min.cluster:
        shutil.move(log_name, 'logs')

    return True


def lmps_step_md(s, bonds, attempt, settings):
    """pysimm.apps.polymatic.lmps_step_md

    Runs LAMMPS step md for the Polymatic algorithm.

    Args:
        s: :class:`~pysimm.system.System` to minimize
        bonds: number of bond to be made
        attempt: number of bonding attempt
        settings: object containing Polymatic settings

    Returns:
        result from :func:`~pysimm.lmps.md`
    """
    if settings.polym.step.cluster:
        nanohub = {'cores': int(settings.polym.step.nanohub_cores),
                   'walltime': int(settings.polym.step.nanohub_walltime)}
        log_name = 'step_%03d_%03d' % (bonds, attempt)
    else:
        nanohub = {}
        log_name = 'logs/step_%03d_%03d' % (bonds, attempt)

    if settings.polym.step.user_input:
        sim = lmps.Simulation(s, name='bond %s attempt #%d' % (bonds + 1, attempt),
                 print_to_screen=False, nanohub=nanohub, custom=True)
        sim.add(settings.polym.step.step_in)
        sim.run(np=settings.np, nanohub=nanohub)
    else:
        sim = lmps.Simulation(s, name='bond %s: attempt #%d' % (bonds + 1, attempt),
            print_to_screen=False, nanohub=nanohub,
            log=log_name
        )
        sim.add(lmps.Init(cutoff=settings.polym.step.nb_cutoff, forcefield=settings.forcefield))
        sim.add(lmps.Velocity(temperature=settings.polym.step.temp))
        sim.add_md(
            ensemble='nvt', temperature=settings.polym.step.temp,
            run=settings.polym.step.length,
        )
        sim.run(np=settings.np, nanohub=nanohub)

    if settings.polym.step.cluster:
        shutil.move(log_name, 'logs')

    return True


def lmps_cycle_nvt_md(s, bonds, settings):
    """pysimm.apps.polymatic.lmps_cycle_nvt_md

    Runs LAMMPS nvt cycle md for the Polymatic algorithm.

    Args:
        s: :class:`~pysimm.system.System` to minimize
        bonds: number of bond to be made
        settings: object containing Polymatic settings

    Returns:
        result from :func:`~pysimm.lmps.md`
    """
    if settings.polym.cycle_nvt.cluster:
        nanohub = {'cores': int(settings.polym.cycle_nvt.nanohub_cores),
                   'walltime': int(settings.polym.cycle_nvt.nanohub_walltime)}
        log_name = 'cycle_nvt_%03d' % bonds
    else:
        nanohub = {}
        log_name = 'logs/cycle_nvt_%03d' % bonds

    if settings.polym.cycle_nvt.user_input:
        sim = lmps.Simulation(s, name='bond %d cycle nvt' % bonds,
                 print_to_screen=False, nanohub=nanohub, custom=True)
        sim.add(settings.polym.cycle_nvt.step_in)
        sim.run(np=settings.np, nanohub=nanohub)
    else:
        sim = lmps.Simulation(s, name='bond %d cycle nvt' % bonds,
            print_to_screen=False, nanohub=nanohub,
            log=log_name
        )
        sim.add(lmps.Init(cutoff=settings.polym.cycle_nvt.nb_cutoff, forcefield=settings.forcefield))
        sim.add(lmps.Velocity(temperature=settings.polym.cycle_nvt.temp))
        sim.add_md(
            ensemble='nvt', temperature=settings.polym.cycle_nvt.temp,
            run=settings.polym.cycle_nvt.length,
        )
        sim.run(np=settings.np, nanohub=nanohub)

    if settings.polym.cycle_nvt.cluster:
        shutil.move(log_name, 'logs')

    return True


def lmps_cycle_npt_md(s, bonds, settings):
    """pysimm.apps.polymatic.lmps_cycle_npt_md

    Runs LAMMPS npt cycle md for the Polymatic algorithm.

    Args:
        s: :class:`~pysimm.system.System` to minimize
        bonds: number of bond to be made
        settings: object containing Polymatic settings

    Returns:
        result from lmps.md
    """
    if settings.polym.cycle_npt.cluster:
        nanohub = {'cores': int(settings.polym.cycle_npt.nanohub_cores),
                   'walltime': int(settings.polym.cycle_npt.nanohub_walltime)}
        log_name = 'cycle_npt_%03d' % bonds
    else:
        nanohub = {}
        log_name = 'logs/cycle_npt_%03d' % bonds

    if settings.polym.cycle_npt.user_input:
        sim = lmps.Simulation(s, name='bond %d cycle npt' % bonds,
                 print_to_screen=False, nanohub=nanohub, custom=True)
        sim.add(settings.polym.cycle_npt.step_in)
        sim.run(np=settings.np, nanohub=nanohub)
    else:
        sim = lmps.Simulation(s, name='bond %d cycle npt' % bonds,
            print_to_screen=False, nanohub=nanohub,
            log=log_name
        )
        sim.add(lmps.Init(cutoff=settings.polym.cycle_npt.nb_cutoff, forcefield=settings.forcefield))
        sim.add(lmps.Velocity(temperature=settings.polym.cycle_npt.temp))
        sim.add_md(
            ensemble='npt', temperature=settings.polym.cycle_npt.nb_cutoff,
            run=settings.polym.cycle_npt.length,
            pressure=settings.polym.cycle_npt.pressure,
        )
        sim.run(np=settings.np, nanohub=nanohub)

    if settings.polym.cycle_npt.cluster:
        shutil.move(log_name, 'logs')

    return True
