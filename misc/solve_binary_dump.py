import sys, petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc as P
import numpy as np

OptDB = P.Options()
fname_mat = OptDB.getString('-mat', 'mat.bin')
fname_b   = OptDB.getString('-b', 'b.bin')

print("Loading Matrix File: {}".format(fname_mat))

viewer = P.Viewer().createBinary(fname_mat, 'r')
A = P.Mat().load(viewer)
A = A.setUp()

viewer = P.Viewer().createBinary(fname_b, 'r')
b = P.Vec().load(viewer)
b = b.setUp()

x = b.duplicate()

ksp = P.KSP().create()
ksp.setOperators(A)
ksp.setFromOptions()
ksp.solve(b, x)
