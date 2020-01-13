# ******************************************************************************
# pysimm.appps.zeoplusplus module
# ******************************************************************************
#
# api to zeoplusplus simulation code
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2019 Ping Lin, Coray M. Colina
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

import sys, os
from subprocess import Popen, PIPE
from time import strftime
import shlex

try:
    from Rappture.tools import getCommandOutput as RapptureExec
except ImportError:
    pass

try:
    from pysimm import system
except ImportError:
    print("Pysimm is required to process pysimm system")
    pass

ZEOpp_EXEC = os.environ.get('ZEOpp_EXEC')

def network(s, **kwargs):
    """pysimm.apps.zeopp.network

    Perform 1. Pore diameters; 2. Channel identification and dimensionality; 3. Surface area;
            4. Accessible volume; 5. Pore size distribution calculation using zeo++ v2.2

    with options to do 6. Probe-occupiable volume; 7. Stochastic ray tracing; 8. Blocking spheres;
                       9. Distance grids; 10. Structure analysis

    Args:
        s: pysimm System object or filename of file in CSSR | CUC | V1 | CIF format
        atype_name: True to use atom type as atom name (usually need radii and mass info), False to use atom element
        radii: file name that contain atom radii data (rad.rad)
        mass: file name that contain atom mass data (mass.mass)
        probe_radius: radius of a probe used in sampling of surface (1.2 A)
        chan_radius: radius of a probe used to determine accessibility of void space (1.2 A)
        num_samples: number of Monte Carlo samples per unit cell (50000)
        option to include in the simulation: set True to activate
            ha: default=True, for using high accuracy,
            res: default=True, for diameters of the largest included sphere, the largest free sphere and the largest included sphere along free sphere path
            chan: default=True, for channel systems characterized by dimensionality as well as Di, Df and Dif
            sa: default=True, for surface area accessible to a spherical probe, characterized by
                      accessible surface area (ASA) and non-accessible surface area (NASA)
            vol: default=True, for accessible volume (AV) and non-accessible volume (NAV)
            volpo: default=False, for accessible proce-occupiable volume (POAV) and non-accessible probe-occupiable volume (PONAV)
            psd: default=True, for the "deriviative distribution" (change of AV w.r.t probe size) reported in the histogram file with 1000 bins of size of 0.1 Ang
            ray_atom: default=False
            block: default=False
            extra: user provided options, such as -gridG, -gridBOV, -strinfo, -oms, etc.
            
    ZEOpp_EXEC: path to zeo++ executable (network)

    Returns:
        None
    """
    global ZEOpp_EXEC

    if ZEOpp_EXEC is None:
        print('Please specify the environment variable ''ZEOpp_EXEC'' that points to '
              'zeo++ executable (network)')
        exit(1)

    probe_radius = kwargs.get('probe_radius', 1.2)
    chan_radius = kwargs.get('chan_radius', 1.2)
    num_samples = kwargs.get('num_samples', 50000)
    atype_name = kwargs.get('atype_name', False)

    ha = kwargs.get('ha', True)
    res = kwargs.get('res', True)
    chan = kwargs.get('chan', True)
    sa = kwargs.get('sa', True)
    vol = kwargs.get('vol', True)
    psd = kwargs.get('psd', True)
    volpo = kwargs.get('volpo', False)
    ray_atom = kwargs.get('ray_atom', False)
    block = kwargs.get('block', False)
    extra = kwargs.get('extra')
 
    nanohub = kwargs.get('nanohub')

    if isinstance(s, system.System):
        if atype_name:
            s.write_cssr('zeopp_data.cssr', aname=1)
        else:
            s.write_cssr('zeopp_data.cssr')
        input_file = 'zeopp_data.cssr'
    elif isinstance(s, str):
        input_file = s
        
    args = ZEOpp_EXEC

    if 'radii' in kwargs.keys(): 
        args += ' -r ' + kwargs.get('radii')
    if 'mass' in kwargs.keys(): 
        args += ' -mass ' + kwargs.get('mass')

    if ha:
        args += ' -ha'
    if res:
        args += ' -res'
    if chan:
        args += ' -chan ' + str(probe_radius) 
    if sa:
        args += ' -sa ' + str(chan_radius) + ' ' + str(probe_radius) + ' ' + str(num_samples) 
    if vol:
        args += ' -vol ' + str(chan_radius) + ' ' + str(probe_radius) + ' ' + str(num_samples) 
    if psd:
        args += ' -psd ' + str(chan_radius) + ' ' + str(probe_radius) + ' ' + str(num_samples) 
    if volpo:
        args += ' -volpo ' + str(chan_radius) + ' ' + str(probe_radius) + ' ' + str(num_samples)
    if ray_atom:
        args += ' -ray_atom ' + str(chan_radius) + ' ' + str(probe_radius) + ' ' + str(num_samples)
    if block:
        args += ' -block ' + str(probe_radius) + ' ' + str(num_samples)
    if extra:
        args += ' ' + extra

    args += ' ' + input_file

    arg_list = shlex.split(args)

    print('%s: starting simulation using zeo++'
          % strftime('%H:%M:%S'))

    if nanohub:
        print('%s: sending zeo++ simulation to computer cluster' % strftime('%H:%M:%S'))
        sys.stdout.flush()
        cmd = ('submit -n 1 -w %s ' % (24*60)) + ZEOpp_EXEC + args
        cmd = shlex.split(cmd)
        exit_status, stdo, stde = RapptureExec(cmd)
    else:
        p = Popen(arg_list, stdin=PIPE, stdout=PIPE, stderr=PIPE) 
        while True:
            stout = p.stdout.readline()
            if stout == '' and p.poll() is not None:
                break
            if stout:
                print(stout.strip())
        # print(stout)
        sterr = p.stderr.readlines()
        print(sterr)

    print('%s: zeo++ simulation successful'
          % strftime('%H:%M:%S'))


