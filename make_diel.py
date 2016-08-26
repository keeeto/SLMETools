#!/usr/bin/env python

import sys
from numpy import array
from numpy import array, diag

with open('EPS', 'r') as f:
    eps = f.readlines()
nedos = len(eps)/6  # this actually equals nedos+1 but it makes things easy

def parse_sarray(data):
    """Helper function for parsing the EPS file.

    Args:
        The lines of the EPS file for a specfic component of the EPS in the form

    Returns:
        The data converted into floats
    """
    return [[float(i) for i in l.split()] for l in data]

xx = array(parse_sarray(eps[1:nedos]))
yy = array(parse_sarray(eps[nedos*1+1:nedos*2]))
zz = array(parse_sarray(eps[nedos*2+1:nedos*3]))
xy = array(parse_sarray(eps[nedos*3+1:nedos*4]))
yz = array(parse_sarray(eps[nedos*4+1:nedos*5]))
zx = array(parse_sarray(eps[nedos*5+1:nedos*6]))
energies = array(xx[:, 0])
real_zip = array(zip(xx[:, 1], yy[:, 1], zz[:, 1],
                     xy[:, 1], yz[:, 1], zx[:, 1]))
imag_zip = array(zip(xx[:, 2], yy[:, 2], zz[:, 2],
                     xy[:, 2], yz[:, 2], zx[:, 2]))

# VASP seems to fail pretty badly at diagonalising the EPS matrix using this
# method, so lets do it ourselves (thanks Alexey Sokol and John Buckeridge)
real_avg = []
imag_avg = []
for energy, real, imag in zip(energies, real_zip, imag_zip):
    real_mat = [[real[0], real[3], real[5]],
                [real[3], real[1], real[4]],
                [real[5], real[4], real[2]]]
    imag_mat = [[imag[0], imag[3], imag[5]],
                [imag[3], imag[1], imag[4]],
                [imag[5], imag[4], imag[2]]]
    real_avg.append(diag(real_mat))
    imag_avg.append(diag(imag_mat))

imag_string = """  frequency dependent IMAGINARY DIELECTRIC FUNCTION (independent particle, no local field effects)
     E(ev)      X         Y         Z        XY        YZ        ZX
       --------------------------------------------------------------------------------------
       """


real_string="""
   frequency dependent      REAL DIELECTRIC FUNCTION (independent particle, no local field effects)
     E(ev)      X         Y         Z        XY        YZ        ZX
  --------------------------------------------------------------------------------------
"""

with open("DIEL", 'w') as f:
    f.write(imag_string)

    for e, imag in zip(energies, imag_avg):
        f.write("%.8f    %.8f    %.8f    %.8f    0.00000      0.00000     0.00000\n" % (e, imag[0], imag[1], imag[2]))

    f.write(real_string)

    for e, real in zip(energies, real_avg):
        f.write("%.8f    %.8f    %.8f    %.8f   0.00000      0.00000     0.00000\n" % (e, real[0], real[1], real[2]))

