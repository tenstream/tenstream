#!/usr/bin/env python

from mpi4py import MPI
from py_boxmc import m_py_boxmc as BMC
import numpy as np

comm = MPI.COMM_WORLD
fcomm = comm.py2f()
myid = comm.Get_rank()

print('This is myid {}'.format(myid))

optprop = [1e-3, 1e-30, .5]
src = 1
ldir = True
phi0 = 0
theta0 = 60
dx, dy, dz = 100, 100, 100

atol, rtol = 1e-2, 1e-1

for theta0 in np.arange(0,60,5):
    for src in np.arange(1,9):
        S, T, Stol, Ttol = BMC.get_coeff_8_10(fcomm, optprop, src, ldir, phi0, theta0, dx, dy, dz, atol, rtol)
        print '{:2d} {:2d}'.format(theta0,src),' '.join(['{:6.3f}'.format(_) for _ in T ]), '::', np.sum(T)

#print(T)
#print(S)

MPI.Finalize()
