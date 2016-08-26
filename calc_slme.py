#!/usr/bin/env python
from operator import itemgetter
from numpy import *
from scipy import linalg
import sys
import re

##############################################################
#
# Reference: L.Yu and A. Zunger, Phys. Rev. Lett. 108, 068701 (2012)
#
##############################################################
#-------------------------------------------------------------------------------------------------------------------
# usage: calc_slme.py sys.argv[1] sys.argv[2] sys.argv[3] sys.argv[4]
# input parameters:
# sys.argv[1]: name of file that contains real and imaginal dielectric permitivity (e.g., copy and paste from OUTCAR)
# sys.argv[2]: minimum band gap, i.e., minEg
# sys.argv[3]: direct allowed band gap, i.e., Egda
# sys.agrv[4]: film thickness in unit of cm
#--------------------------------------------------------------------------------------------------------------------
diel = open(sys.argv[1],'r').readlines()
Eg = float(sys.argv[2])                     # minEg
dEg =  float(sys.argv[3])  - Eg             # Eg^{da}-minEg
L  = float(sys.argv[4])     # thin film thickness (in unit of cm): e.g, 2.0E-4 cm

am15 = open("AM15-eta_SQ.dat",'r').readlines() # solar spectrum
ev = 1.60217648740E-19
h = 6.626068E-34
c = 299792458
Pin = 1000/ev

#
# calculating absorption efficient
#
f=open('alpha_w.dat','w')
f.write("#  E(hv)   alpha_x       alpha_y       alpha_z       alpha_av    n_av**2\n")

npts = (len(diel) - 8)/2
sigma = 0.06/2.3548
dE = float(diel[5].split()[0]) - float(diel[4].split()[0])

alpha = []

eps_imag = []
lthres = 0

for npt in range(npts) :
    xi  = diel[npt+4].split()
    xi1 = diel[npt+5].split()
    eps_imag.append([float(xi[1]),float(xi[2]),float(xi[3]), float(xi[4]), float(xi[5]),float(xi[6])])
    if float(xi[0]) <= 5.0  and float(xi1[0]) >= 5.0 :
       nmax = npt
#    if float(xi[1])+float(xi[2])+float(xi[3]) > 0.0 and lthres == 0:
#       thres_Eg = float(xi[0])
#       lthres = 1


for npt in range(nmax) :

    xr = diel[npt+8+npts].split()
    hv = float(xr[0])
    R_eps_xx = float(xr[1])
    R_eps_yy = float(xr[2])
    R_eps_zz = float(xr[3])
    R_eps_xy = float(xr[4])
    R_eps_yz = float(xr[5])
    R_eps_zx = float(xr[6])


    I_eps_xx = 0.0
    I_eps_yy = 0.0
    I_eps_zz = 0.0
    I_eps_xy = 0.0
    I_eps_yz = 0.0
    I_eps_zx = 0.0
    # gaussian broadening  eps_imag
    for npt2 in range(npts) :
        val = ((npt-npt2)*dE/sigma)**2/2.0
        if val < 15.0 :
           wt = exp(-val)
           I_eps_xx = I_eps_xx + wt*eps_imag[npt2][0]
           I_eps_yy = I_eps_yy + wt*eps_imag[npt2][1]
           I_eps_zz = I_eps_zz + wt*eps_imag[npt2][2]
           I_eps_xy = I_eps_xy + wt*eps_imag[npt2][3]
           I_eps_yz = I_eps_yz + wt*eps_imag[npt2][4]
           I_eps_zx = I_eps_zx + wt*eps_imag[npt2][5]

    I_eps_xx = I_eps_xx*dE/sqrt(2.0*pi)/sigma
    I_eps_yy = I_eps_yy*dE/sqrt(2.0*pi)/sigma
    I_eps_zz = I_eps_zz*dE/sqrt(2.0*pi)/sigma
    I_eps_xy = I_eps_xy*dE/sqrt(2.0*pi)/sigma
    I_eps_yz = I_eps_yz*dE/sqrt(2.0*pi)/sigma
    I_eps_zx = I_eps_zx*dE/sqrt(2.0*pi)/sigma

    Cxx = complex(R_eps_xx, I_eps_xx)
    Cyy = complex(R_eps_yy, I_eps_yy)
    Czz = complex(R_eps_zz, I_eps_zz)
    Cxy = complex(R_eps_xy, I_eps_xy)
    Cyz = complex(R_eps_yz, I_eps_yz)
    Czx = complex(R_eps_zx, I_eps_zx)

    C_eps = mat([[Cxx, Cxy, conj(Czx)],[conj(Cxy), Cyy, Cyz], [Czx, conj(Cyz), Czz]])

    eps_eig, eps_v = linalg.eig(C_eps)

#   print hv, imag(eps_eig[0]),imag(eps_eig[1]),imag(eps_eig[2])
    alpha_a1 = hv * 71618.96076 * sqrt(abs(eps_eig[0])-real(eps_eig[0]))
    alpha_a2 = hv * 71618.96076 * sqrt(abs(eps_eig[1])-real(eps_eig[1]))
    alpha_a3 = hv * 71618.96076 * sqrt(abs(eps_eig[2])-real(eps_eig[2]))
    alpha_av = (alpha_a1 + alpha_a2 + alpha_a3)/3

    n1 = sqrt(0.5*(abs(eps_eig[0])+real(eps_eig[0])))
    n2 = sqrt(0.5*(abs(eps_eig[1])+real(eps_eig[1])))
    n3 = sqrt(0.5*(abs(eps_eig[2])+real(eps_eig[2])))
    n_av = (n1+n2+n3)/3.0

    f.write("%9.4f %13.6E %13.6E %13.6E %13.6E %9.4f\n" %(hv, alpha_a1, alpha_a2, alpha_a3, alpha_av, n_av**2))
    alpha.append([hv, alpha_av, n_av**2])

#
# preparing the data for calculating slme
#
data_in = []

for l in range(1, len(am15)) :
    x = am15[l].split()
    hv  = float(x[0])
    nhv = float(x[1])

    for ll in range(len(alpha)-1) :
        if alpha[ll][0] <= hv and alpha[ll+1][0] >= hv :
           fact = (hv - alpha[ll][0])/(alpha[ll+1][0] - alpha[ll][0])
           tmp1 = alpha[ll][1]*(1-fact) + fact*alpha[ll+1][1]
           tmp2 = alpha[ll][2]*(1-fact) + fact*alpha[ll+1][2]
           data_in.append([hv, nhv, tmp1, tmp2])
           break


# calculate Isc and I0

Vc = 0.025851997434
Isc = 0.0
I0  = 0.0

for l in range(len(data_in)-1) :
    hv0 = data_in[l][0]
    hv1 = data_in[l+1][0]
    #
    des = hv1 - hv0
    #
    aE0 = 1.0 - exp(-2.0*L*data_in[l][2])
    aE1 = 1.0 - exp(-2.0*L*data_in[l+1][2])

    is0 = data_in[l][1]*aE0
    is1 = data_in[l+1][1]*aE1

    Isc = Isc + (is0 + is1)*des/2.0

    irb0 = hv0**2/(exp(hv0/Vc)-1)*aE0
    irb1 = hv1**2/(exp(hv1/Vc)-1)*aE1
#   irb0 = hv0**2/(exp(hv0/Vc)-1)*aE0*data_in[l][3]
#   irb1 = hv1**2/(exp(hv1/Vc)-1)*aE1*data_in[l+1][3]

    I0 = I0 + (irb0 + irb1)*des/2.0

I0 = I0 * 2.0*pi*c/(h*c/ev)**3 * exp(dEg/Vc)

#
# calculate max IV = [ Isc - I0*e^{(V+dEg)/Vc} ] *V
# dE = 0.001
#
#npts = int(thres_Eg/0.001)
npts = int(Eg/0.001)

maxIV = 0
for ll in range(npts) :
    Vap = ll*0.001
    IVtmp = Vap * ( Isc - I0*exp(Vap/Vc))
    if IVtmp > maxIV :
       maxIV = IVtmp
       #%Vm = Vap

slme = maxIV/Pin*100.0

print slme

f.close()
