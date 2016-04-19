# ******************************************************************************
# pysimm.apps.glass_transition module
# ******************************************************************************
#
# tools to extract tg
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

import numpy as np
from itertools import izip
import numpy as np
from matplotlib import pyplot as plt
from pysimm import system, lmps, calc
from pysimm import error_print

rappture = True
try:
    import Rappture
except ImportError:
    rappture = False


def chunks(l, n):
    # Yield successive n-sized chunks from l.
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def read_log_file(file_='cooling.log'):
    step = []
    temp = []
    density = []
    with file(file_) as f:
        for line in f:
            if line.strip() and line.split()[0] == 'Step':
                line = f.next().split()
                while line[0] != 'Loop':
                    step.append(int(line[0]))
                    temp.append(float(line[1]))
                    density.append(float(line[2]))
                    line = f.next().split()

    return step, temp, density


def log_avg(step, temp, density, ncool, nequil, t_range=None):
    avg_temp = []
    std_temp = []
    avg_den = []
    std_den = []
    itemp = iter(temp)
    iden = iter(density)
    for t in t_range:
        ttemp = []
        tden = []
        for n in range(nequil+1):
            ttemp.append(itemp.next())
            tden.append(iden.next())
        avg_temp.append(np.mean(ttemp))
        avg_den.append(np.mean(tden))
        std_temp.append(np.std(ttemp))
        std_den.append(np.std(tden))
        if t == t_range[-1]:
            break
        else:
            for n in range(ncool+1):
                itemp.next()
                iden.next()
    return avg_temp, avg_den, std_temp, std_den


def write_den_data(step, temp, density, outfile='density.dat'):
    with file(outfile, 'w') as f:
        for s, t, d in izip(step, temp, density):
            f.write('{}\t{}\t{}\n'.format(s, t, d))


def polyfit(x, y, degree):
    # solution from http://stackoverflow.com/questions/893657/how-do-i-calculate-r-squared-using-python-and-numpy

    coeffs = np.polyfit(x, y, degree)

    p = np.poly1d(coeffs)
    yhat = p(x)
    ybar = np.sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((y - ybar)**2)
    return coeffs, ssreg / sstot


def find_fitting_range(temp, den, plot=False):
    nlower = 0
    best_r2 = None
    best_coeffs = None
    for x in range(3, len(temp)):
        coeffs, r2 = polyfit(temp[:x], den[:x], 1)
        p = np.poly1d(coeffs)
        y = p(temp[:x])
        if y[-1] >= den[:x][-1]:
            nlower += 1
        else:
            nlower = 0
        if nlower < 5:
            best_coeffs, best_r2 = coeffs, r2
        if plot:
            plt.plot(temp[:x], den[:x], temp[:x], y)
            plt.show()
    return best_coeffs, best_r2


def fit_and_plot(temp, den, temperr, denerr):
    hi_coeffs, hi_r2 = find_fitting_range(temp, den)
    lo_coeffs, lo_r2 = find_fitting_range(list(reversed(temp)), list(reversed(den)))
    lo_fit = np.poly1d(lo_coeffs)
    hi_fit = np.poly1d(hi_coeffs)
    plt.errorbar(temp, den, xerr=temperr, yerr=denerr)
    plt.plot(temp, lo_fit(temp), temp, hi_fit(temp))
    print calc.intersection(((0, hi_fit(0)), (1, hi_fit(1))), ((0, lo_fit(0)), (1, lo_fit(1))))
    plt.show()
    return calc.intersection(((0, hi_fit(0)), (1, hi_fit(1))), ((0, lo_fit(0)), (1, lo_fit(1))))


def calc_tg(logfile='cooling.log', ncool=None, nequil=None, t_range=None):
    if not ncool or not nequil or not t_range:
        error_print('must provide number of cooling data points (ncool), number of equilibration data points(nequil) '
                    'and temperature range(t_range)')
    fit_and_plot(*log_avg(*read_log_file(logfile), ncool=ncool, nequil=nequil, t_range=t_range))


def cool(s, **kwargs):
    hi_t = kwargs.get('hi_t') or 1000
    lo_t = kwargs.get('lo_t') or 100
    delta_t = kwargs.get('delta_t') or 100
    cooling_length = kwargs.get('cooling_length') or 600000
    equil_length = kwargs.get('equil_length') or 100000

    thermo = kwargs.get('thermo') or 1000
    print_to_screen = kwargs.get('print_to_screen')

    np = kwargs.get('np')

    sim = lmps.Simulation(s, name='cooling', print_to_screen=print_to_screen)
    sim.add_md(lmps.MolecularDynamics(temp=hi_t, ensemble='npt', pressure=1, length=equil_length, new_v=True,
                                      thermo=thermo, thermo_style='custom step temp density'))
    for t in reversed(range(lo_t+delta_t, hi_t+1, delta_t)):
        sim.add_md(lmps.MolecularDynamics(t_start=t, t_stop=t-delta_t, ensemble='npt', pressure=1,
                                          length=cooling_length, scale_v=True, thermo=thermo,
                                          thermo_style='custom step temp density'))
        sim.add_md(lmps.MolecularDynamics(temp=t-delta_t, ensemble='npt', pressure=1,
                                          length=equil_length, scale_v=True, thermo=thermo,
                                          thermo_style='custom step temp density'))

    sim.run(np=np)
