# ******************************************************************************
# pysimm.apps.equilibrate module
# ******************************************************************************
#
# 21-step equilibration algorithm written using pysimm tools
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2016 Michael E. Fortunato
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from itertools import izip

from pysimm import system, lmps

rappture = True
try:
    import Rappture
except ImportError:
    rappture = False


def equil(s, **kwargs):
    """pysimm.apps.equilibrate.equil

    Runs a 21-step compression/decompression equilibration algorithm

    Args:
        s: :class:`~pysimm.system.System` object
        tmax: maximum temperature during equilibration
        pmax: maximum pressure during equilibration
        tfinal: desired final temperature of final system
        pfinal: desired final pressure of final system
        np: number of processors to use during equilibration simulations
        p_steps: list of pressures to use during equilibration (must match length of length_list)
        length_list: list of simulation durations to use during equilibration (must match length of p_steps)

    Returns:
        None
    """
    tmax = kwargs.get('tmax', 1000)
    pmax = kwargs.get('pmax', 50000)
    tfinal = kwargs.get('tfinal', 300)
    pfinal = kwargs.get('pfinal', 1)
    
    init = kwargs.get('init')
    output_settings = kwargs.get('output_settings')

    np = kwargs.get('np')

    p_list = kwargs.get('p_steps', [0.02*pmax, 0.6*pmax, pmax, 0.5*pmax, 0.1*pmax, 0.01*pmax, pfinal])

    length_list = kwargs.get('length_list', [100000, 100000, 100000, 100000, 100000, 100000, 100000])

    sim = lmps.Simulation(s, name='equil', **kwargs)
    
    if init:
        sim.add(init)
    if output_settings:
        sim.add(output_settings)
        
    sim.add(lmps.Velocity(temperature=tfinal))

    step = 0
    for p, l in izip(p_list, length_list):
        step += 1
        if l:
            sim.add_md(length=l/2, ensemble='nvt', temperature=tmax, **kwargs)
            sim.add_md(length=l, ensemble='nvt', temperature=tfinal, **kwargs)
            sim.add_md(length=l/2, ensemble='npt', temperature=tfinal, pressure=p, **kwargs)

    sim.run(np=np)

    s.write_lammps('equil.lmps')
    s.write_xyz('equil.xyz')
