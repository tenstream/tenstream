#!/usr/bin/env python
""" Create Input for TenStream caluclations on a gaussian hill"""

import numpy as np

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


def create_hill_x_input():
    """ Create Input for TenStream caluclations on a gaussian hill"""
    from scipy.ndimage.filters import gaussian_filter as gf
    from scipy.ndimage.filters import gaussian_filter1d as gf1d
    from scipy.interpolate import interp1d

    # Create elevation map for gaussian hill
    Nx, Ny = 3, 31
    
    DX = .1  # 100m horizontal resolution
    HILL_HEIGHT = 1.5  # 500m hill

    elev = np.zeros((Nx, Ny))
    elev[:,(Ny-1)/2] = 1.
    elev = gf1d(elev, HILL_HEIGHT/DX, axis=1)
    elev /= np.max(elev) / HILL_HEIGHT

    # Interpolate atmosphere file on new grid
    afglus = np.loadtxt('afglus_100m.dat')

    z = afglus[:,0]
    p = afglus[:,1]

    Nz = np.shape(z)[-1]

    pressure_by_height = interp1d(z, p, kind='linear')

    hhl = np.tile(z, Nx*Ny).reshape((Nx, Ny, Nz))[:, :, ::-1]  # hhl now begins at surface 
    hill = hhl.copy()

    # Create surface following sigma coordinates with coordinate stretching in the vertical
    for k in range(Nz):
        hill[:, :, k] = hhl[:, :, k] + elev * np.exp(-hhl[:, :, k] / HILL_HEIGHT)

    hill = np.minimum(hill, np.max(z))
    hill = hill[:,:,::-1]  # hhl now begins at TOA
    p_hill = pressure_by_height(hill)

    lev_coord = p_hill

    create_var('plev'  , interp_var(p, afglus[:, 1], lev_coord), 'lev')

    lay_coord = (lev_coord[:,:,1:] + lev_coord[:,:,:-1]) / 2

    create_var('tlay',   interp_var(p, afglus[:, 2], lay_coord), 'lay')
    create_var('air',    interp_var(p, afglus[:, 3], lay_coord), 'lay')
    create_var('o3vmr',  interp_var(p, afglus[:, 4], lay_coord), 'lay')
    create_var('o2vmr',  interp_var(p, afglus[:, 5], lay_coord), 'lay')
    create_var('h2ovmr', interp_var(p, afglus[:, 6], lay_coord), 'lay')
    create_var('co2vmr', interp_var(p, afglus[:, 7], lay_coord), 'lay')
    create_var('n2ovmr', interp_var(p, afglus[:, 8], lay_coord), 'lay')

if __name__ == '__main__':
    create_hill_x_input()
