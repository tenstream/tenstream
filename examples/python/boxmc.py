#!/usr/bin/env python

from mpi4py import MPI
from py_boxmc import m_py_boxmc as BMC

comm = MPI.COMM_WORLD
fcomm = comm.py2f()
myid = comm.Get_rank()

print('This is myid {}'.format(myid))

optprop = [1e-3, 1e-3, .5]
src = 1
ldir = True
phi0 = 0
theta0 = 0
dx, dy, dz = 100, 100, 50

atol, rtol = 1e-2, 1e-1

S, T, Stol, Ttol = BMC.get_coeff_8_10(fcomm, optprop, src, ldir, phi0, theta0, dx, dy, dz, atol, rtol)

print(T)
print(S)
