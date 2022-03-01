#!/usr/bin/env python
import sys, petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc as P
import numpy as np
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
myid = comm.Get_rank()

OptDB = P.Options()
fname_mat = OptDB.getString('-mat', 'mat.bin')
fname_b   = OptDB.getString('-b', 'b.bin')

check_matrows = OptDB.getBool('-check_matrows', False)
info = OptDB.getBool('-info', False)

if myid==0:
    print("Loading Matrix File: {}".format(fname_mat))

viewer = P.Viewer().createBinary(fname_mat, 'r')
A = P.Mat().load(viewer)
A = A.setUp()

if check_matrows:
    if myid==0:
        print("Building Transpose...")
    At = A.duplicate()
    A.transpose(At)
    if myid==0:
        print("Building Transpose... done")

    rStart, rEnd = At.getOwnershipRange()
    for i in range(rStart, rEnd):
        col_idx, coeff = At.getRow(i)
        if np.sum(coeff) < -1e-8:
            print("Have column vector with sum lt 1: icol", i, np.where(coeff!=0), \
                    ':', coeff, 'sum', np.sum(coeff))
        else:
            if info:
                print("Col Vec {:8d} seems OK:".format(i), ' '.join(['{:+0.2f}'.format(c) for c in coeff]))


if check_matrows:
    if myid==0:
        print("Checking Matrix rows if they are OK")
    rStart, rEnd = A.getOwnershipRange()
    for i in range(rStart, rEnd):
        col_idx, coeff = A.getRow(i)
        if not np.isfinite(coeff).all():
            print("Have NaN in mat coeffs", coeff)

        if any(coeff>1.):
            print("Have coeffs gt 1", coeff)

        if any(coeff<-1):
            print("Have coeffs lt -1", coeff)

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

if myid==0:
    print("Solving System")
solve_start = time.time()
ksp.solve(b, x)
if myid==0:
    print('Solve took about {} s'.format(time.time()-solve_start))

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
