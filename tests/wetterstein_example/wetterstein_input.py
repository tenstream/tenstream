#!/usr/bin/env python
""" Create Input for TenStream//rrtmg caluclations surrounding the Wetterstein massif """

import numpy as np


def load_srtm_wetterstein(xmin=3000, xmax=3200, ymin=1100, ymax=1450):
    """ Load elevation data from SRTM surrounding the Wetterstein massif """
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open('./srtm_39_03.tif')
    elev = ds.ReadAsArray()
    return (elev[xmin:xmax, ymin:ymax]/1e3).T


def create_var(varname, axis, vals, newaxis, dim_prefix):
    import os
    import netCDF4 as NC
    from scipy.interpolate import interp1d

    fp = interp1d(axis, vals, kind='linear', fill_value='extrapolate')
    var = fp(newaxis)[:, :, ::-1]

    fname = 'wetterstein_input.nc'
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

def create_srtm_input():
    """ Create Input for TenStream//rrtmg caluclations surrounding the Wetterstein massif """
    from scipy.ndimage.filters import gaussian_filter as gf
    from scipy.interpolate import interp1d

    wetterstein = load_srtm_wetterstein(xmin=3000, xmax=3200, ymin=1100, ymax=1300)
    wetterstein = gf(wetterstein,2, mode='wrap')

    Nx, Ny = np.shape(wetterstein)

    afglus = np.loadtxt('afglus_100m.dat')
    z = afglus[:,0]
    Nz = np.shape(z)[-1]

    p = afglus[:,1]

    fp = interp1d(z, p, kind='linear', fill_value='extrapolate')

    hhl = np.tile(z, Nx*Ny).reshape((Nx, Ny, Nz))[:, :, ::-1]  # hhl now begins at surface 
    hill = hhl.copy()

    # Create surface following sigma coordinates with coordinate stretching in the vertical
    for k in range(Nz):
        hill[:, :, k] = hhl[:, :, k] + wetterstein * np.exp(-hhl[:, :, k] / 4.)

    create_var('plev'  , afglus[:, 0], afglus[:, 1], hill, 'lev')

    create_var('tlay',   afglus[:, 0], afglus[:, 2], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('air',    afglus[:, 0], afglus[:, 3], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('o3vmr',  afglus[:, 0], afglus[:, 4], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('o2vmr',  afglus[:, 0], afglus[:, 5], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('h2ovmr', afglus[:, 0], afglus[:, 6], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('co2vmr', afglus[:, 0], afglus[:, 7], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')
    create_var('n2ovmr', afglus[:, 0], afglus[:, 8], (hill[:, :, 1:]+hill[:, :, :-1])*.5, 'lay')


if __name__ == '__main__':
    create_srtm_input()
