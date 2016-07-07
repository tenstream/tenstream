#!/usr/bin/env python
""" Create Input for TenStream caluclations on a gaussian hill"""

import numpy as np


def create_var(varname, axis, vals, newaxis, dim_prefix):
    """ Write regridded data to netcdf file """
    import os
    import netCDF4 as NC
    from scipy.interpolate import interp1d

    fp = interp1d(axis, vals, kind='linear', fill_value = "extrapolate" )
    var = fp(newaxis)[:, :, ::-1]

    fname = 'hill1_x_input.nc'
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

def create_hill_x_input():
    """ Create Input for TenStream caluclations on a gaussian hill"""
    from scipy.ndimage.filters import gaussian_filter as gf
    from scipy.ndimage.filters import gaussian_filter1d as gf1d
    from scipy.interpolate import interp1d

    # Create elevation map for gaussian hill
    Nx, Ny = 3, 16
    
    elev = np.zeros((Nx, Ny))
    elev[:,(Ny-1)/2] = 1.
    elev = gf1d(elev, (Ny-1)/4, axis=1) 

    # Interpolate atmosphere file on new grid
    afglus = np.loadtxt('afglus_100m.dat')
    z = afglus[:,0]
    Nz = np.shape(z)[-1]

    p = afglus[:,1]

    fp = interp1d(z, p, kind='linear')

    hhl = np.tile(z, Nx*Ny).reshape((Nx, Ny, Nz))[:, :, ::-1]  # hhl now begins at surface 
    hill = hhl.copy()

    # Create surface following sigma coordinates with coordinate stretching in the vertical
    for k in range(Nz):
        hill[:, :, k] = hhl[:, :, k] + elev * np.exp(-hhl[:, :, k] / 4.)

    create_var('plev'  , afglus[:, 0], afglus[:, 1], hill, 'lev')

    create_var('tlay',   afglus[:, 0], afglus[:, 2], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('air',    afglus[:, 0], afglus[:, 3], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('o3vmr',  afglus[:, 0], afglus[:, 4], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('o2vmr',  afglus[:, 0], afglus[:, 5], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('h2ovmr', afglus[:, 0], afglus[:, 6], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('co2vmr', afglus[:, 0], afglus[:, 7], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('n2ovmr', afglus[:, 0], afglus[:, 8], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')


if __name__ == '__main__':
    create_hill_x_input()
