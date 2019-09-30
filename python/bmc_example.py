#!/bin/env python3
from pylab import *
from mpi4py import MPI
from python.py_boxmc import m_py_boxmc as P

def delta_scale(kabs, ksca, g, f=None, max_g=None):
    import numpy as np
    if f is None:
        f = g**2

    if max_g is not None:
        f = (max_g-g) / (max_g-1)
        if g < max_g:
            return np.array([kabs, ksca, g])

    dtau = kabs+ksca
    w0   = ksca/dtau

    dtau_d = dtau * ( 1. - w0 * f )
    g_d    = ( g - f ) / ( 1. - f )
    w0_d   = w0 * ( 1. - f ) / ( 1. - f * w0 )

    kabs_d = dtau_d * (1. - w0_d)
    ksca_d = dtau_d * w0_d

    return np.array([kabs_d, ksca_d, g_d])

comm = MPI.COMM_WORLD.py2f()
dx=100.
dz=100.
kabs=1e-10
ksca=1e-1
g=0.
src=1
direct=True
phi = 0
theta = 0
atol=1e-3
rtol=5e-1
S, T, S_tol, T_tol = P.py_get_coeff_8_10(comm, dx, dx, dz, [kabs, ksca, g], src, direct, phi, theta, atol, rtol)

print("S",S)
print("T",T)
print("S_tol",S_tol)
print("T_tol",T_tol)


# Test wedge coeffs
g_range = linspace(0,.85,30)

ret            = [ P.py_get_coeff_wedge_5_8(comm, dx, dx, dz,            [kabs, ksca, g]           , src, direct, phi, theta, atol, rtol) for g in g_range ]
ret_delta      = [ P.py_get_coeff_wedge_5_8(comm, dx, dx, dz, delta_scale(kabs, ksca, g)           , src, direct, phi, theta, atol, rtol) for g in g_range ]
ret_delta_maxg = [ P.py_get_coeff_wedge_5_8(comm, dx, dx, dz, delta_scale(kabs, ksca, g, max_g=.65), src, direct, phi, theta, atol, rtol) for g in g_range ]

def unpack_results(x):
    S, T, S_tol, T_tol = [ np.array([ g_entry[i] for g_entry in x ]) for i in range(4) ]
    S[:,-1] += np.sum(T, axis=-1)
    return S, T, S_tol, T_tol

#plt.figure(1)
#plt.clf()
#S, T, Stol, Ttol = unpack_results(ret)
#
#[ plt.plot(g_range, S[:,i], label="S dst={}".format(i)) for i in range(len(S[0])) ]
#xlabel('assym. param g')
#plt.legend(loc='best')
#
#plt.figure(2)
#plt.clf()
#S, T, Stol, Ttol = unpack_results(ret_delta)
#
#[ plt.plot(g_range, S[:,i], label="delta scaled S dst={}".format(i)) for i in range(len(S[0])) ]
#xlabel('assym. param g')
#plt.legend(loc='best')
#
#plt.figure(3)
#plt.clf()
#S, T, Stol, Ttol = unpack_results(ret_delta_maxg)
#
#[ plt.plot(g_range, S[:,i], label="delta scaled S dst={}".format(i)) for i in range(len(S[0])) ]
#xlabel('assym. param g')
#plt.legend(loc='best')

plt.figure(num=10, figsize=(20,12))
plt.clf()
for i in range(len(ret[0][0])):
    sub=plt.subplot(2,4,i+1)

    S, T, Stol, Ttol = unpack_results(ret)
    plt.plot(g_range, S[:,i], label="S dst={}".format(i))

    S, T, Stol, Ttol = unpack_results(ret_delta)
    plt.plot(g_range, S[:,i], label="delta S dst={}".format(i))

    S, T, Stol, Ttol = unpack_results(ret_delta_maxg)
    plt.plot(g_range, S[:,i], label="delta_maxg S dst={}".format(i))

    plt.legend(loc='best')
plt.suptitle('BoxMC Wedge 5_8')
plt.tight_layout()
plt.savefig('delta_scaling_experiment_bmc_wedge_5_8_dir2diff_src{}.pdf'.format(src), bbox_inches='tight')


run_bmc_8_10 = lambda op_func: [ np.mean([P.py_get_coeff_8_10(comm, dx, dx, dz, op_func([kabs, ksca, g]), isrc, direct, phi, theta, atol, rtol) for isrc in range(1,5) ], axis=0) for g in g_range ]

ret8_10            = run_bmc_8_10( lambda x:x )
ret8_10_delta      = run_bmc_8_10( lambda x: delta_scale(*x) )
ret8_10_delta_maxg = run_bmc_8_10( lambda x: delta_scale(*x, max_g=.65) )


def unpack_results(x):
    S, T, S_tol, T_tol = [ np.array([ g_entry[i] for g_entry in x ]) for i in range(4) ]
    S[:,1] += np.sum(T, axis=-1)
    return S, T, S_tol, T_tol

plt.figure(num=11, figsize=(20,12))
plt.clf()
for i in range(len(ret8_10[0][0])):
    sub=plt.subplot(3,4,i+1)

    S, T, Stol, Ttol = unpack_results(ret8_10)
    plt.plot(g_range, S[:,i], label="S dst={}".format(i))

    S, T, Stol, Ttol = unpack_results(ret8_10_delta)
    plt.plot(g_range, S[:,i], label="delta S dst={}".format(i))

    S, T, Stol, Ttol = unpack_results(ret8_10_delta_maxg)
    plt.plot(g_range, S[:,i], label="delta_maxg S dst={}".format(i))

    plt.legend(loc='best')
plt.suptitle('BoxMC 8_10')
plt.tight_layout()
plt.savefig('delta_scaling_experiment_bmc_8_10_dir2diff_src{}.pdf'.format(src), bbox_inches='tight')



run_bmc_3_10 = lambda op_func: [ P.py_get_coeff_3_10(comm, dx, dx, dz, op_func([kabs, ksca, g]), src, direct, phi, theta, atol, rtol) for g in g_range ]

ret3_10            = run_bmc_3_10( lambda x:x )
ret3_10_delta      = run_bmc_3_10( lambda x: delta_scale(*x) )
ret3_10_delta_maxg = run_bmc_3_10( lambda x: delta_scale(*x, max_g=.65) )


def unpack_results(x):
    S, T, S_tol, T_tol = [ np.array([ g_entry[i] for g_entry in x ]) for i in range(4) ]
    S[:,1] += np.sum(T, axis=-1)
    return S, T, S_tol, T_tol

plt.figure(num=12, figsize=(20,12))
plt.clf()
for i in range(len(ret3_10[0][0])):
    sub=plt.subplot(3,4,i+1)

    S, T, Stol, Ttol = unpack_results(ret3_10)
    plt.plot(g_range, S[:,i], label="S dst={}".format(i))

    S, T, Stol, Ttol = unpack_results(ret3_10_delta)
    plt.plot(g_range, S[:,i], label="delta S dst={}".format(i))

    S, T, Stol, Ttol = unpack_results(ret3_10_delta_maxg)
    plt.plot(g_range, S[:,i], label="delta_maxg S dst={}".format(i))

    plt.legend(loc='best')
plt.suptitle('BoxMC 3_10')
plt.tight_layout()
plt.savefig('delta_scaling_experiment_bmc_3_10_dir2diff_src{}.pdf'.format(src), bbox_inches='tight')


direct=False
src=2
run_bmc_10 = lambda op_func: [ P.py_get_coeff_wedge_5_8(comm, dx, dx, dz, op_func([kabs, ksca, g]), src, direct, phi, theta, atol, rtol) for g in g_range ]

ret10            = run_bmc_10( lambda x:x )
ret10_delta      = run_bmc_10( lambda x: delta_scale(*x) )
ret10_delta_maxg = run_bmc_10( lambda x: delta_scale(*x, max_g=.65) )


def unpack_results(x):
    S, T, S_tol, T_tol = [ np.array([ g_entry[i] for g_entry in x ]) for i in range(4) ]
    return S, T, S_tol, T_tol

plt.figure(num=13, figsize=(20,12))
plt.clf()
for i in range(len(ret10[0][0])):
    sub=plt.subplot(2,4,i+1)

    S, T, Stol, Ttol = unpack_results(ret10)
    plt.plot(g_range, S[:,i], label="S dst={}".format(i))

    S, T, Stol, Ttol = unpack_results(ret10_delta)
    plt.plot(g_range, S[:,i], label="delta S dst={}".format(i))

    S, T, Stol, Ttol = unpack_results(ret10_delta_maxg)
    plt.plot(g_range, S[:,i], label="delta_maxg S dst={}".format(i))

    plt.legend(loc='best')
plt.suptitle('BoxMC Wedge 5_8 Diffuse src={}'.format(src))
plt.tight_layout()
plt.savefig('delta_scaling_experiment_bmc_wedge_5_8_diff2diff_src{}.pdf'.format(src), bbox_inches='tight')
