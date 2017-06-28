#!/usr/bin/env python

import numpy as np
from scipy.interpolate import interp1d

from mpi4py import MPI
import py_rrtm_lw_sw as RRTM

TS = RRTM.m_py_rrtm_lw_sw


def load_cld(fname):
    """
        Function to read a libRadtran Cloud file
    """
    from numpy import loadtxt,arange,zeros,array
    print('opening cloud file',fname)
    f = open(fname,'r')
    Nx,Ny,Nz = array(f.readline().split()).astype(int)[:3]
    hhl = array(f.readline().split()).astype(float)
    f.close()
    geometry = {
            'Nx' : Nx,
            'Ny' : Ny,
            'Nz' : Nz,
            'dx' : hhl[0]*1e3,
            'dy' : hhl[1]*1e3,
            'z'  : hhl[2:]*1e3,
            }
    x,y,z,v1,v2 = loadtxt(fname,skiprows=2,unpack=True)
    x=x.astype(int)
    y=y.astype(int)
    z=z.astype(int)
    lwc  = zeros((Nz,Ny,Nx ))
    reff = zeros((Nz,Ny,Nx ))
    for i in arange(len(v1)):
            lwc  [z[i]-1,y[i]-1,x[i]-1] = v1[i]
            reff [z[i]-1,y[i]-1,x[i]-1] = v2[i]
    return array([lwc, reff]), geometry


def py_rrtmg(
        comm=MPI.COMM_WORLD,
        Nx=3, Ny=3, Nz=10, zt=None, max_height=1e4,
        dx=100, dy=100, phi0=0, theta0=0,
        albedo_thermal=0, albedo_solar=.15,
        atm_filename='afglus_100m.dat',
        lthermal=True, lsolar=True,
        N_ranks_y=None, N_ranks_x=None,
        lwc=None, reliq=None, iwc=None, reice=None):

    fcomm = comm.py2f()
    myid = comm.Get_rank()
    numnodes = comm.Get_size()

    if (N_ranks_x is None) or (N_ranks_y is None):
        N_ranks_y = int(np.sqrt(1.*numnodes))
        N_ranks_x = numnodes // N_ranks_y
        if N_ranks_y*N_ranks_x != numnodes:
            N_ranks_x = numnodes
            N_ranks_y = 1

    nxproc = np.ones(N_ranks_x) * Nx
    nyproc = np.ones(N_ranks_y) * Ny

    solar_albedo_2d = np.ones((Nx, Ny)) * albedo_solar
    opt_time = 0
    icollapse = 1

    if zt is None:
        zt = np.linspace(0, max_height, Nz+1)

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

    if lwc is None:
        lwc   = np.ones_like(co2) * 0
    if reliq is None:
        reliq = np.ones_like(co2) * 4
    if iwc is None:
        iwc   = np.ones_like(co2) * 0
    if reice is None:
        reice = np.ones_like(co2) * 20

    # Either call with reduced interface
    TS.rrtmg_minimal(fcomm,
            dx, dy, phi0, theta0,
            albedo_thermal, albedo_solar,
            atm_filename,
            lthermal, lsolar,
            plev, tlev,
            lwc, reliq, iwc, reice,
            nxproc, nyproc)

    edir, edn, eup, abso = [ x.copy() for x in [TS.edir, TS.edn, TS.eup, TS.abso]]

    TS.destroy_rrtmg()

    def stitch(local):
        """ Gather values From local to global grid """
        g = comm.gather(local, root=0)
        if myid == 0:
            return np.concatenate([np.concatenate(g[i*N_ranks_x:(i+1)*N_ranks_x],axis=1) for i in range(N_ranks_y) ], axis=2)
        return None

    edir, edn, eup, abso = [ stitch(x) for x in [edir, edn, eup, abso]]

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


def small_example(theta=60, phi=180, hhl=3e3, Nz=20, cld_z=15):
    """
        Small example with a liquid water cloud in at along the x-axis.
        Size of a rank's subdomain is constant(i.e. more nodes -> bigger domain)

    """
    myid = MPI.COMM_WORLD.Get_rank()
    nproc = MPI.COMM_WORLD.Get_size()
    dx = dy = 5e2
    Nx,Ny,Nz = 3, 15, Nz
    lwc, reff = np.zeros((Nz,Nx,Ny)), np.ones((Nz,Nx,Ny))*10
    lwc[cld_z, :, Ny//4] = 1

    edir,edn,eup,abso = py_rrtmg(Nx=lwc.shape[1], Ny=lwc.shape[2], Nz=lwc.shape[0],
            dx=dx, dy=dy, max_height=hhl,
            lwc=lwc,
            theta0=theta, phi0=phi,
            N_ranks_x=1, N_ranks_y=MPI.COMM_WORLD.Get_size(),
            lthermal=False)

    import matplotlib.pyplot as plt
    plt.figure(1); plt.clf()
    plt.subplot(211)
    plt.imshow(edir[-Nz:,Nx//2,:],interpolation='nearest', extent=(0,Ny*dy,0,hhl))
    plt.colorbar()
    plt.subplot(212)
    plt.plot(np.linspace(0,Ny*dy, Ny), edir[-cld_z,Nx//2,:], label=r"$E_{dir} @ cloudlayer$")
    plt.plot(np.linspace(0,Ny*dy, Ny), edir[-1,Nx//2,:], label=r"$E_{dir} @ srfc$")
    plt.plot(np.linspace(0,Ny*dy, Ny), edn [-1,Nx//2,:], label=r"$E_{dn}  @ srfc$")
    plt.plot(np.linspace(0,Ny*dy, Ny), eup [-1,Nx//2,:], label=r"$E_{up}  @ srfc$")
    plt.legend(framealpha=.5)
    plt.waitforbuttonpress()

    return edir,edn,eup,abso


def i3rc_example(theta=60, phi=180):
    """
        Computes the radiative transfer for a cloud field from the radiative transfer intercomparison project.
        Domain decomposition happens only in y-direction and number of processes has to be divisor of 96
    """
    myid = MPI.COMM_WORLD.Get_rank()
    nproc = MPI.COMM_WORLD.Get_size()

    (lwc, reff), geom = load_cld('wcloud.i3rc1.dat')
    Nx, Ny, Nz, dx, dy, hhl = [ geom[k] for k in ('Nx', 'Ny', 'Nz', 'dx', 'dy', 'z') ]

    ymin, ymax = Ny * np.array([myid, myid+1])/nproc
    myslice = slice(ymin, ymax, None)

    lwc = lwc[:, :,  myslice]
    reff = reff[:, :, myslice]

    edir,edn,eup,abso = py_rrtmg(Nx=lwc.shape[1], Ny=lwc.shape[2], Nz=lwc.shape[0], max_height=hhl[-1],
            lwc=lwc,
            theta0=theta, phi0=phi,
            N_ranks_x=1, N_ranks_y=MPI.COMM_WORLD.Get_size(),
            lthermal=False)

    if myid == 0:
        np.save('i3rc1_{}'.format(theta), (edir,edn,eup,abso))
        import matplotlib.pyplot as plt
        plt.clf()
        plt.subplot(141); plt.title('edir')
        plt.imshow(edir[:,0],vmin=0)
        plt.colorbar()
        plt.subplot(142); plt.title('edn')
        plt.imshow(edn[:,0],vmin=0)
        plt.colorbar()
        plt.subplot(143); plt.title('eup')
        plt.imshow(eup[:,0],vmin=0)
        plt.colorbar()
        plt.subplot(144); plt.title('abso')
        plt.imshow(abso[:,0]*86)
        plt.colorbar()

        [sub.label_outer() for sub in plt.gcf().get_axes()]
        plt.gcf().set_size_inches((7.4,8))
        plt.tight_layout()

        plt.savefig('i3rc1_{}.pdf'.format(theta))
        return edir,edn,eup,abso


if __name__ == "__main__":
    small_example()
    pass
