#!/usr/bin/env python

import numpy as np
from scipy.interpolate import interp1d

from mpi4py import MPI
from py_rrtm_lw_sw import m_py_rrtm_lw_sw as RRTM

def py_rrtmg(
        comm=MPI.COMM_WORLD,
        Nx=3, Ny=3, Nz=10,
        dx=100, dy=100,
        phi0=0, theta0=0,
        albedo_thermal=0,
        albedo_solar=.15,
        lthermal=True,
        lsolar=True):

    fcomm = comm.py2f()
    myid = comm.Get_rank()
    numnodes = comm.Get_size()

    N_ranks_y = int(np.sqrt(1.*numnodes))
    N_ranks_x = numnodes // N_ranks_y
    if N_ranks_y*N_ranks_x != numnodes:
        N_ranks_x = numnodes
        N_ranks_y = 1

    nxproc = np.ones(N_ranks_x) * Nx
    nyproc = np.ones(N_ranks_y) * Ny

    if myid == 0:
        print(myid, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y, '::', numnodes)

    solar_albedo_2d = np.ones((Nx, Ny)) * albedo_solar
    opt_time = 0
    icollapse = 1

    zt = np.linspace(0, 1e4, Nz+1)

    atm_filename = 'afglus_100m.dat'

    z, p, T, air, o3, o2, h2o, co2, no2 = np.loadtxt(atm_filename, unpack=True)
    z *= 1e3
    o3, o2, h2o, co2, no2 = [ x/air for x in [o3, o2, h2o, co2, no2]]

    afglus2dynamics = lambda new_x, old_x, var: np.repeat( interp1d(old_x, var)(new_x) ,Nx*Ny).reshape((len(new_x), Nx,Ny))

    plev = afglus2dynamics(zt, z, p)
    play = (plev[1:] + plev[:-1]) / 2

    iplev = plev[:,0,0]
    iplay = play[:,0,0]

    tlev = afglus2dynamics(iplev, p, T)

    tlay, o3, o2, h2o, co2, no2 = [afglus2dynamics(iplay, p, x) for x in [T, o3, o2, h2o, co2, no2]]

    ch4   = np.ones_like(co2) * 2e-6
    lwc   = np.ones_like(co2) * 0
    reliq = np.ones_like(co2) * 4
    iwc   = np.ones_like(co2) * 0
    reice = np.ones_like(co2) * 20

    # Either call with reduced interface
    RRTM.rrtmg_minimal(fcomm,
            dx, dy, phi0, theta0,
            albedo_thermal, albedo_solar,
            atm_filename,
            lthermal, lsolar,
            plev, tlev,
            lwc, reliq, iwc, reice)

    edir, edn, eup, abso = [ x.copy() for x in [RRTM.edir, RRTM.edn, RRTM.eup, RRTM.abso]]

    RRTM.destroy_rrtmg()

    return edir, edn, eup, abso

# Or provide the full set of options
#RRTM.rrtmg(fcomm,
#        dx, dy, phi0, theta0,
#        albedo_thermal, albedo_solar,
#        atm_filename,
#        lthermal, lsolar,
#        plev, tlev, nxproc, nyproc,
#        tlay, h2o, o3,
#        co2, ch4, no2,  o2,
#        lwc, reliq, iwc, reice,
#        icollapse, opt_time,
#        solar_albedo_2d)
#
## Results are then available as RRTM module attributes
#if myid == 0:
#    for k, (dr, dn, up) in enumerate(zip(RRTM.edir, RRTM.edn, RRTM.eup)):
#        print('{:4d} {:12.2f} {:12.2f} {:12.2f}'.format(k, dr[0,0], dn[0,0], up[0,0]))


#while True:
#   edir, edn, eup, abso = py_rrtmg()


#MPI.Finalize()
