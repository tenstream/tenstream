#!/usr/bin/env python
import sys, petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc as P
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
myid = comm.Get_rank()

OptDB = P.Options()
fname_mat = OptDB.getString('-mat', 'mat.bin')
fname_b   = OptDB.getString('-b', 'b.bin')

check_matrows = OptDB.getBool('-check_matrows', False)

if myid==0:
    print("Loading Matrix File: {}".format(fname_mat))

viewer = P.Viewer().createBinary(fname_mat, 'r')
A = P.Mat().load(viewer)
A = A.setUp()
if check_matrows:
    rStart, rEnd = A.getOwnershipRange()
    for i in range(rStart, rEnd):
        col_idx, coeff = A.getRow(i)
        if not np.isfinite(coeff).all():
            print("Have NaN in mat coeffs", coeff)

        if any(coeff>1.):
            print("Have coeffs gt 1", coeff)

        if any(coeff<-1):
            print("Have coeffs lt -1", coeff)

        cv = A.getColumnVector(i)
        if np.sum(cv.array)<0:
            print("Have column vector with sum lt 1: icol", i, np.where(cv.array!=0), \
                    ':', cv.array[cv.array!=0]), 'sum', np.sum(cv.array)

if myid==0:
    print("Loading src Vec: {}".format(fname_b))
viewer = P.Viewer().createBinary(fname_b, 'r')
b = P.Vec().load(viewer)
b = b.setUp()
if not np.isfinite(b.array).all():
    print(myid,"Have NaN in src term")

x = b.duplicate()

ksp = P.KSP().create()
ksp.setOperators(A)
ksp.setFromOptions()
ksp.solve(b, x)

fname_plot = OptDB.getString('-plot', '')
if fname_plot is not '':
    from pylab import *
    Ncells = 27
    Nlevel = 179
    solution_top = x.array[0:Ncells*Nlevel:Nlevel]
    solution_bot = x.array[Nlevel-1:Ncells*Nlevel:Nlevel]
    subplot(211); plot(solution_top)
    subplot(212); plot(solution_bot)
    savefig(fname_plot)
