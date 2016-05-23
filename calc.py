# ******************************************************************************
# pysimm.calc module
# ******************************************************************************
#
# vector rotation
# particle separation distance
# particle separation distance considering PBC
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

from random import random
from itertools import izip
from math import sin, cos, pi, acos
import numpy as np

from pysimm import error_print
from pysimm import warning_print
from pysimm import verbose_print
from pysimm import debug_print
from pysimm.utils import Item
from pysimm.utils import ItemContainer


def hyperbola(x, a, b, c, d, e):
    #   hyperbola(x) with parameters
    #   a/b = asymptotic slope
    #    c  = curvature at vertex
    #    d  = offset to vertex
    #    e  = vertical offset
    return a*np.sqrt((b*c)**2 + (x-d)**2)/b + e


def rot_hyperbola(x, a, b, c, d, e, th):
    pars = a, b, c, 0, 0  # do the shifting after rotation
    xd = x - d
    hsin = hyperbola(xd, *pars)*np.sin(th)
    xcos = xd*np.cos(th)
    return e + hyperbola(xcos - hsin, *pars)*np.cos(th) + xcos - hsin


def intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y


def find_rotation(a, b):
    a = np.array(a)
    b = np.array(b)

    a_x_b = np.cross(a, b)
    axis = a_x_b / np.linalg.norm(a_x_b)

    theta = acos(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))

    skew = np.matrix([[0, -axis[2], axis[1]],
                     [axis[2], 0, -axis[0]],
                     [-axis[1], axis[0], 0]])

    rot_matrix = np.identity(3) + sin(theta)*skew + (1 - cos(theta))*skew*skew

    return rot_matrix


def rotate_vector(x, y, z, theta_x=None, theta_y=None, theta_z=None):
    xt = random() * 2 * pi if theta_x is None else theta_x
    yt = random() * 2 * pi if theta_y is None else theta_y
    zt = random() * 2 * pi if theta_z is None else theta_z

    c = np.matrix([[x], [y], [z]])

    rot_mat_x = np.matrix([[1, 0, 0],
                           [0, cos(xt), -sin(xt)],
                           [0, sin(xt), cos(xt)]])
    rot_mat_y = np.matrix([[cos(yt), 0, sin(yt)],
                           [0, 1, 0],
                           [-sin(yt), 0, cos(yt)]])
    rot_mat_z = np.matrix([[cos(zt), -sin(zt), 0],
                           [sin(zt), cos(zt), 0],
                           [0, 0, 1]])

    c = rot_mat_x * c
    c = rot_mat_y * c
    c = rot_mat_z * c

    return [x[0] for x in c.tolist()]


def distance(p1, p2):
    return np.linalg.norm([p1.x - p2.x, p1.y - p2.y, p1.z - p2.z])


def angle(p1, p2, p3, radians=False):
    p12 = distance(p1, p2)
    p23 = distance(p2, p3)
    p13 = distance(p1, p3)
    theta = acos((pow(p12, 2)+pow(p23, 2)-pow(p13, 2))/(2*p13*p23))
    if not radians:
        theta = theta * 180 / pi
    return theta


def chiral_angle(a, b, c, d):
    ht = np.array([a.x-b.x, a.y-b.y, a.z-b.z])
    ht /= np.linalg.norm(ht)

    hmethyl = np.array([a.x-d.x, a.y-d.y, a.z-d.z])
    hmethyl /= np.linalg.norm(hmethyl)

    hside = np.array([a.x-c.x, a.y-c.y, a.z-c.z])
    hside /= np.linalg.norm(hside)

    side_x_methyl = np.cross(hside, hmethyl)

    side_dot_methyl = np.dot(hside, hmethyl)
    side_theta_methyl = acos(side_dot_methyl)

    side_x_methyl /= sin(side_theta_methyl)

    cos_theta = np.dot(side_x_methyl, ht)

    return acos(cos_theta)/pi*180


def auto_tacticity(s, return_angles=True, unwrap=True, rewrap=True):

    # buggy, because chirality has to do with chains.................

    s.add_particle_bonding()

    stereochem_angles = ItemContainer()

    if unwrap:
        s.unwrap()

    for p in s.particles:
        if p.chiral is True:
            bonded_mw = {}
            for pb in p.bonded_to:
                bonded_mw[pb] = pb.type.mass
                for pb_ in pb.bonded_to:
                    if pb_ is not p:
                        bonded_mw[pb] += pb_.type.mass
            print [(p_.tag, bonded_mw[p_]) for p_ in bonded_mw]
            print [x.tag for x in sorted(bonded_mw, key=bonded_mw.get)]
            sorted_by_mw = sorted(bonded_mw, key=bonded_mw.get)
            stereochem_angles.add(Item(tag=p.tag, value=chiral_angle(p,
                                                                     sorted_by_mw[0],
                                                                     sorted_by_mw[1],
                                                                     sorted_by_mw[2])))

    last = None
    iso_diads = 0
    syn_diads = 0
    for a in stereochem_angles:
        if last is not None:
            if (a.value < 90 and last.value < 90) or (a.value > 90 and last.value > 90):
                iso_diads += 1
            else:
                syn_diads += 1
        last = a

    if iso_diads == (len(stereochem_angles) - 1):
        t = 'isotactic'
    elif syn_diads == (len(stereochem_angles) - 1):
        t = 'syndiotactic'
    else:
        t = 'atactic'

    if rewrap:
        s.wrap()

    if return_angles:
        return t, stereochem_angles
    else:
        return t


def tacticity(s, a_tag=None, b_tag=None, c_tag=None, d_tag=None, offset=None, return_angles=True, unwrap=True,
              rewrap=True, skip_first=False):

    if a_tag is None or b_tag is None or c_tag is None or d_tag is None:
        error_print('particle tags for chiral center are required')
        error_print('a: chiral center, b-d: 3 side groups')

    if offset is None:
        error_print('offset for tags in each monomer is required, i.e. - number of particles in each monomer')

    if unwrap:
        s.unwrap()

    a_ = [s.particles[i] for i in range(a_tag, s.particles.count, offset)]
    b_ = [s.particles[i] for i in range(b_tag, s.particles.count, offset)]
    c_ = [s.particles[i] for i in range(c_tag, s.particles.count, offset)]
    d_ = [s.particles[i] for i in range(d_tag, s.particles.count, offset)]

    stereochem_angles = []

    if skip_first:
        for a, b, c, d in izip(a_[1:], b_[1:], c_[1:], d_[1:]):
            stereochem_angles.append(chiral_angle(a, b, c, d))
    else:
        for a, b, c, d in izip(a_, b_, c_, d_):
            stereochem_angles.append(chiral_angle(a, b, c, d))

    last = None
    iso_diads = 0
    syn_diads = 0
    for a in stereochem_angles:
        if last is not None:
            if (a < 90 and last < 90) or (a > 90 and last > 90):
                iso_diads += 1
            else:
                syn_diads += 1
        last = a

    if iso_diads == (len(stereochem_angles) - 1):
        t = 'isotactic'
    elif syn_diads == (len(stereochem_angles) - 1):
        t = 'syndiotactic'
    else:
        t = 'atactic'

    if rewrap:
        s.wrap()

    if return_angles:
        return t, stereochem_angles
    else:
        return t


def frac_free_volume(v_sp, v_void):
    return (-0.3 * v_sp + 1.3 * v_void)/v_sp


def pbc_distance(s, p1, p2):
    frac_x1 = p1.x / s.dim.dx
    frac_y1 = p1.y / s.dim.dy
    frac_z1 = p1.z / s.dim.dz

    frac_x2 = p2.x / s.dim.dx
    frac_y2 = p2.y / s.dim.dy
    frac_z2 = p2.z / s.dim.dz

    frac_d = np.array([frac_x1 - frac_x2, frac_y1 - frac_y2, frac_z1 - frac_z2])
    frac_d = frac_d - np.round(frac_d)

    dx = frac_d[0] * s.dim.dx
    dy = frac_d[1] * s.dim.dy
    dz = frac_d[2] * s.dim.dz

    return np.linalg.norm([dx, dy, dz])
