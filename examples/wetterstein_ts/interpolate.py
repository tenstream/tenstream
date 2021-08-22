# The script interpolates *_rf.nc files hhl_rf, abso_rf and lwc_rf from the sigma coordinate mesh onto a regular mesh for easier plotting in visit

import xarray as xr
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import griddata
from scipy.interpolate import Rbf
import os
import argparse

parser = argparse.ArgumentParser(description='Interpolate *_rf.nc hhl_rf, abso_rf and lwc_rf from sigma coordinate mesh onto regular coordinates and save as *_rf_itpd.nc')
parser.add_argument('-f', '--fnames-list', nargs='+', default=[])
fnames = parser.parse_args().fnames_list

print('List of files to interpolate:', fnames)

for fname in fnames:

    print('Starting with {}'.format(fname))

    if (os.path.exists(fname)):
        f = xr.open_dataset(fname)

        z = f['hhl_rf'].data[:]
        z = ((z[:,:,:-1] + z[:,:,1:]) / 2) * 1e-3
        xx = np.tile(np.arange(int(f.Nx/2)), (int(f.Ny/2), 1)) * 100 * 1e-3
        yy = np.tile(np.arange(int(f.Ny/2)), (int(f.Nx/2), 1)).T * 100 * 1e-3
        abso = f['abso_rf'].data[:]
        lwc = f['lwc_rf'].data[:]

        x = np.empty(z.shape)
        y = np.empty(z.shape)
        for i in range(z.shape[2]):
            x[:,:,i] = xx
            y[:,:,i] = yy

        c = np.empty([z.flatten().shape[0],3])
        c[:,0] = x.flatten()
        c[:,1] = y.flatten()
        c[:,2] = z.flatten()

#rbf =  Rbf(x,y,z,abso)

        xtrgt = np.arange(int(f.Nx/2)) * 100 * 1e-3
        ytrgt = np.arange(int(f.Ny/2)) * 100 * 1e-3
        ztrgt = (np.arange(100) / 99 * np.max(z))[::-1]
        X, Y, Z = np.meshgrid(xtrgt, ytrgt, ztrgt)

        ctrgt = np.empty([X.flatten().shape[0],3])
        ctrgt[:,0] = X.flatten()
        ctrgt[:,1] = Y.flatten()
        ctrgt[:,2] = Z.flatten()

        dfLin = griddata(c, abso.flatten(), ctrgt, method='nearest')
        #dfNear = griddata(c, abso.flatten(), ctrgt, method='nearest')
        #dfLin[np.isnan(dfLin)] = 0#dfNear[np.isnan(dfLin)]

        abso_out = dfLin.reshape(X.shape)

        dfLin = griddata(c, lwc.flatten(), ctrgt, method='nearest')
        #dfNear = griddata(c, abso.flatten(), ctrgt, method='nearest')
        #dfLin[np.isnan(dfLin)] = 0#dfNear[np.isnan(dfLin)]

        lwc_out = dfLin.reshape(X.shape)


        D = xr.Dataset({
            'abso': xr.DataArray(abso_out, dims=('dim1','dim2','dim3')),
            'lwc': xr.DataArray(lwc_out, dims=('d3','d2','d1'))
            })
        D.to_netcdf(fname.replace('.nc', '_itpd.nc'))

        print('Done with {}'.format(fname))

    else:
        print('Could not find {}'.format(fname))

