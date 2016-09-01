#!/bin/python
from __future__ import print_function
def spin_polarised(filename='OUTCAR'):
    '''
    Checks whether the calculation is spin polarised or not
    '''
    spin = False
    for line in open(filename):
            if len(line.split()) > 5 and line.split()[0] == 'ISPIN':
                if float(line.split()[2]) == 2:
                    spin = True
    return spin 

def number_of_electrons(filename='OUTCAR'):
    nelect = None
    for line in open('OUTCAR'):
        if line.find('total number of electrons') != -1:
            nelect = float(line.split('=')[1].split()[0].strip())
    return nelect

def number_of_bands(filename='OUTCAR'):
    bands = None
    for line in open('OUTCAR'):
        if line.find('number of bands') != -1:
            bands = float(line.split()[14])
    return bands

def get_bands(nbands, filename='OUTCAR', spin=1):
    occ = []
    unocc = []
    lines = open(filename).readlines()
    for n, line in enumerate(lines):
        if len(line.split()) == 6 and line.split()[0] == "k-point":
            inp = line.split()
            kpt = [float(inp[3]), float(inp[4]), float(inp[5])]

            for i in range(n+2, n+2+nbands):
                inp = lines[i].split()
                if float(inp[2]) == 0.0:
                    unocc.append([kpt, float(inp[1])])
                else:
                    occ.append([kpt, float(inp[1])])
    return occ, unocc

nbands = int(number_of_bands())
occupied, unoccupied = get_bands(nbands)
occupied = sorted(occupied,key=lambda x: x[1],reverse=True)
unoccupied = sorted(unoccupied,key=lambda x: x[1],reverse=False)

print("VBM: {p[1]:8.2f}. k-point: {o[0]:6.4f} {o[1]:6.4f} {o[2]:6.4f}.".format(p = occupied[0], o = occupied[0][0]))
print("CBM: {p[1]:8.2f}. k-point: {o[0]:6.4f} {o[1]:6.4f} {o[2]:6.4f}.".format(p = unoccupied[0], o = unoccupied[0][0]))
