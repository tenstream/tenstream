#!/usr/bin/env python
""" Create Input for TenStream//rrtmg caluclations surrounding the Wetterstein massif """

import numpy as np


def load_srtm_wetterstein():
    """ Load elevation data from SRTM surrounding the Wetterstein massif """
    import netCDF4 as NC
    D = NC.Dataset('means_wetterstein.cdf')
    elev = D['elevation'][255:511, 255:511]
    D.close()
    return elev.T

def create_var(varname, var, dim_prefix):
    """ Write regridded data to netcdf file """
    import os
    import netCDF4 as NC
    from scipy.interpolate import interp1d

    fname = 'input.nc'
    if os.path.exists(fname):
        D = NC.Dataset(fname, 'a')
    else:
        D = NC.Dataset(fname, 'w')


    Nx, Ny, Nz = np.shape(var)

    try:
        D.createDimension(dim_prefix+'Nx', Nx)
        D.createDimension(dim_prefix+'Ny', Ny)
        D.createDimension(dim_prefix+'Nz', Nz)

        x = D.createVariable(dim_prefix+'Nx', 'd', (dim_prefix+'Nx',))
        y = D.createVariable(dim_prefix+'Ny', 'd', (dim_prefix+'Ny',))
        x[:] = np.linspace(0, Nx*.1, Nx)
        y[:] = np.linspace(0, Ny*.1, Ny)
    except:
        pass

    data = D.createVariable(varname, 'd', (dim_prefix+'Ny', dim_prefix+'Nx', dim_prefix+'Nz'))
    data[:] = np.swapaxes(var, 0, 1)
    data[:,:,:] = data[:,:,::-1]

    D.close()


def interp_var(old_pressure_grid, var, new_pressure_grid):
   """ use afglus values and interpolate it on the hill pressure grid """
   from scipy.interpolate import interp1d
   p_us = old_pressure_grid

   Nx, Ny, Nz = np.shape(new_pressure_grid)
   new_var = np.zeros_like(new_pressure_grid) * np.NaN
   for i in range(Nx):
     for j in range(Ny):

       fp = interp1d(p_us, var, kind='linear')

       new_var[i,j] = fp(new_pressure_grid[i,j])
   return new_var


def create_srtm_input():
    """ Create Input for TenStream//rrtmg caluclations surrounding the Wetterstein massif """
    from scipy.ndimage.filters import gaussian_filter as gf
    from scipy.interpolate import interp1d

    # wetterstein = load_srtm_wetterstein(xmin=3000, xmax=3200, ymin=1100, ymax=1300)
    wetterstein = load_srtm_wetterstein()

    # make sure that height values around the border are cyclic
    wetterstein[:,0] = wetterstein[:,-1] = wetterstein[:,[0,-1]].mean(axis=1)
    wetterstein[0,:] = wetterstein[-1,:] = wetterstein[[0,-1],:].mean(axis=0)

    # smooth topography
    wetterstein = gf(wetterstein, 1.5, mode='wrap')

    Nx, Ny = np.shape(wetterstein)

    afglus = np.loadtxt('afglus_100m.dat')
    z = afglus[:,0]
    Nz = np.shape(z)[-1]

    p = afglus[:,1]

    pressure_by_height = interp1d(z, p, kind='linear')

    hhl = np.tile(z, Nx*Ny).reshape((Nx, Ny, Nz))[:, :, ::-1]  # hhl now begins at surface
    hill = hhl.copy()

    # Create surface following sigma coordinates with coordinate stretching in the vertical
    for k in range(Nz):
        hill[:, :, k] = hhl[:, :, k] + wetterstein * np.exp(-hhl[:, :, k] / np.max(wetterstein))

    hill = np.minimum(hill, np.max(z))
    hill = hill[:, :, ::-1]  # hhl now begins at TOA
    p_hill = pressure_by_height(hill)

    kmax = np.argmax(p_hill > 500) # drop everything above 100hPa
    lev_coord = p_hill[:,:,kmax:]

    create_var('plev', interp_var(p, afglus[:, 1], lev_coord), 'lev')
    create_var('tlev', interp_var(p, afglus[:, 2], lev_coord), 'lev')

    if False:
        lay_coord = (lev_coord[:, :, 1:] + lev_coord[:, :, :-1]) / 2

        create_var('tlay',   interp_var(p, afglus[:, 2], lay_coord), 'lay')
        create_var('air',    interp_var(p, afglus[:, 3], lay_coord), 'lay')
        create_var('o3vmr',  interp_var(p, afglus[:, 4], lay_coord), 'lay')
        create_var('o2vmr',  interp_var(p, afglus[:, 5], lay_coord), 'lay')
        create_var('h2ovmr', interp_var(p, afglus[:, 6], lay_coord), 'lay')
        create_var('co2vmr', interp_var(p, afglus[:, 7], lay_coord), 'lay')
        create_var('n2ovmr', interp_var(p, afglus[:, 8], lay_coord), 'lay')


if __name__ == '__main__':
    create_srtm_input()
