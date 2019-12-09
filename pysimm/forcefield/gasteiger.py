# ******************************************************************************
# pysimm.gasteiger module
# ******************************************************************************
#
# gasteiger algorithm and parameters (to be moved to pysimm.forcefield)
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2016 Michael E. Fortunato, Coray M. Colina
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

from .. import error_print
from .. import warning_print
from .. import verbose_print
from .. import debug_print
from ..utils import Item, ItemContainer

element_names_by_mass = {1: 'H', 4: 'He', 7: 'Li', 9: 'Be', 11: 'B', 12: 'C',
                         14: 'N', 16: 'O', 19: 'F', 20: 'Ne', 23: 'Na',
                         24: 'Mg', 27: 'Al', 28: 'Si', 31: 'P', 32: 'S',
                         35: 'Cl', 39: 'K', 40: 'Ca', 80: 'Br', 127: 'I'}

gasteiger_parameters = ItemContainer()
gasteiger_parameters.add(Item(tag='H_', name='H_', a=7.17, b=6.24, c=-0.56))
gasteiger_parameters.add(Item(tag='C_3', name='C_3', a=7.98, b=9.18, c=1.88))
gasteiger_parameters.add(Item(tag='C_2', name='C_2', a=8.79, b=9.32, c=1.51))
gasteiger_parameters.add(Item(tag='C_R', name='C_R', a=8.79, b=9.32, c=1.51))
gasteiger_parameters.add(Item(tag='C_1', name='C_1', a=10.39, b=9.45, c=0.73))
gasteiger_parameters.add(Item(tag='N_3', name='N_3', a=11.54, b=10.82, c=1.36))
gasteiger_parameters.add(Item(tag='N_2', name='N_2', a=12.87, b=11.15, c=0.85))
gasteiger_parameters.add(Item(tag='N_R', name='N_R', a=12.87, b=11.15, c=0.85))
gasteiger_parameters.add(Item(tag='N_1', name='N_1', a=15.68, b=11.7, c=-0.27))
gasteiger_parameters.add(Item(tag='O_3', name='O_3', a=14.18, b=12.92, c=1.39))
gasteiger_parameters.add(Item(tag='O_2', name='O_2', a=17.07, b=13.79, c=0.47))
gasteiger_parameters.add(Item(tag='O_R', name='O_R', a=17.07, b=13.79, c=0.47))
gasteiger_parameters.add(Item(tag='F_', name='F_', a=14.66, b=13.85, c=2.31))
gasteiger_parameters.add(Item(tag='Cl_', name='Cl_', a=11.00, b=9.69, c=1.35))
gasteiger_parameters.add(Item(tag='Br_', name='Br_', a=10.08, b=8.47, c=1.16))
gasteiger_parameters.add(Item(tag='I_', name='I_', a=9.90, b=7.96, c=0.96))
gasteiger_parameters.add(Item(tag='S_', name='S_', a=10.14, b=9.13, c=1.38))
gasteiger_parameters.add(Item(tag='S_3', name='S_3', a=10.14, b=9.13, c=1.38))

gasteiger_parameters.add(Item(tag='h', name='h', a=7.17, b=6.24, c=-0.56))
gasteiger_parameters.add(Item(tag='c_sp3', name='c_sp3', a=7.98, b=9.18,
                              c=1.88))
gasteiger_parameters.add(Item(tag='c_sp2', name='c_sp2', a=8.79, b=9.32,
                              c=1.51))
gasteiger_parameters.add(Item(tag='c_sp1', name='c_sp1', a=10.39, b=9.45,
                              c=0.73))
gasteiger_parameters.add(Item(tag='n_sp3', name='n_sp3', a=11.54, b=10.82,
                              c=1.36))
gasteiger_parameters.add(Item(tag='n_sp2', name='n_sp2', a=12.87, b=11.15,
                              c=0.85))
gasteiger_parameters.add(Item(tag='n_sp1', name='n_sp1', a=15.68, b=11.7,
                              c=-0.27))
gasteiger_parameters.add(Item(tag='o_sp3', name='o_sp3', a=14.18, b=12.92,
                              c=1.39))
gasteiger_parameters.add(Item(tag='o_sp2', name='o_sp2', a=17.07, b=13.79,
                              c=0.47))
gasteiger_parameters.add(Item(tag='f', name='f', a=14.66, b=13.85, c=2.31))
gasteiger_parameters.add(Item(tag='cl', name='cl', a=11.00, b=9.69, c=1.35))
gasteiger_parameters.add(Item(tag='br', name='br', a=10.08, b=8.47,
                              c=1.16))
gasteiger_parameters.add(Item(tag='i', name='i', a=9.90, b=7.96, c=0.96))
gasteiger_parameters.add(Item(tag='s', name='s', a=10.14, b=9.13, c=1.38))


def set_charges(s, maxiter=100, tol=1e-6):
    global gasteiger_parameters

    for p in s.particles:

        if (p.type and p.type.name and p.type.name.find('@') > 0 and
                s.particle_types.get(p.type.name.split('@')[-1])):
            p.nbonds = len(p.bonds) + 1
            if p.type.name[0] == 'H':
                p.linker = 'head'
            elif p.type.name[0] == 'T':
                p.linker = 'tail'
            else:
                p.linker = True
            type_ = s.particle_types.get(p.type.name.split('@')[-1])
            if type_:
                p.type = type_[0]
            else:
                error_print('found linker type %s but did not find regular '
                            'type %s' % (p.type.name,
                                         p.type.name.split('@')[-1]))
                return
        else:
            p.nbonds = len(p.bonds)

        gast_type = gasteiger_parameters.get(p.type.name)
        if gast_type:
            gast_type = gast_type[0]
            p.gast_conv = False
            p.qn = 0.0
            p.charge = 0.0
            p.gast_a = gast_type.a
            p.chi = p.gast_a
            p.gast_b = gast_type.b
            p.gast_c = gast_type.c
            continue
        else:
            if not p.type.elem and p.type.mass:
                elem = element_names_by_mass.get(int(round(p.type.mass)))
                if elem:
                    p.type.elem = elem
            if not p.type.elem or p.type.elem not in ['H', 'N', 'C', 'O',
                                                      'F', 'Cl', 'Br', 'I',
                                                      'S']:
                error_print('cannot find gastieger paramater for particle %s'
                            % p.tag)
            else:
                if p.type.elem == 'H':
                    gast_type = gasteiger_parameters.get('h')[0]
                elif p.type.elem == 'C':
                    if p.nbonds == 4:
                        gast_type = gasteiger_parameters.get('c_sp3')[0]
                    elif p.nbonds == 3:
                        gast_type = gasteiger_parameters.get('c_sp2')[0]
                    elif p.nbonds == 2:
                        gast_type = gasteiger_parameters.get('c_sp1')[0]
                elif p.type.elem == 'N':
                    if p.nbonds >= 3:
                        gast_type = gasteiger_parameters.get('n_sp3')[0]
                    elif p.nbonds == 2:
                        gast_type = gasteiger_parameters.get('c_sp2')[0]
                    elif p.nbonds == 1:
                        gast_type = gasteiger_parameters.get('c_sp1')[0]
                elif p.type.elem == 'O':
                    if p.nbonds == 2:
                        gast_type = gasteiger_parameters.get('o_sp3')[0]
                    elif p.nbonds == 1:
                        gast_type = gasteiger_parameters.get('o_sp2')[0]
                elif p.type.elem == 'F':
                    gast_type = gasteiger_parameters.get('f')[0]
                elif p.type.elem == 'Cl':
                    gast_type = gasteiger_parameters.get('cl')[0]
                elif p.type.elem == 'Br':
                    gast_type = gasteiger_parameters.get('br')[0]
                elif p.type.elem == 'I':
                    gast_type = gasteiger_parameters.get('i')[0]
                elif p.type.elem == 'S':
                    gast_type = gasteiger_parameters.get('s')[0]

                gast_type = gast_type
                p.gast_conv = False
                p.qn = 0.0
                p.charge = 0.0
                p.gast_a = gast_type.a
                p.chi = p.gast_a
                p.gast_b = gast_type.b
                p.gast_c = gast_type.c
                continue

    for n in range(1, maxiter+1):
        for b in s.bonds:
            if b.b.chi > b.a.chi:
                b.b.qn += ((b.a.chi - b.b.chi) /
                           (b.a.gast_a + b.a.gast_b + b.a.gast_c))
                b.a.qn += ((b.b.chi - b.a.chi) /
                           (b.a.gast_a + b.a.gast_b + b.a.gast_c))
            else:
                b.b.qn += ((b.a.chi - b.b.chi) /
                           (b.b.gast_a + b.b.gast_b + b.b.gast_c))
                b.a.qn += ((b.b.chi - b.a.chi) /
                           (b.b.gast_a + b.b.gast_b + b.b.gast_c))

        convergence = True
        for p in s.particles:
            p.qn *= pow(0.5, n)
            p.charge += p.qn
            p.chi = p.gast_a + p.gast_b*p.charge + p.gast_c*p.charge*p.charge
            if abs(p.qn) < tol:
                p.gast_conv = True
            else:
                convergence = False
            p.qn = 0.0

        if convergence:
            print('charges converged after %s iterations' % n)
            break

    if not convergence:
        print('charges not converged after %s iterations' % n)
