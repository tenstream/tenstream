# refolds the wetterstein.F90 output files from quadrupled domain back into single domain

import numpy as np
import xarray as xr
import os
import argparse

parser = argparse.ArgumentParser(description='Interpolate *_rf.nc hhl_rf, abso_rf and lwc_rf from sigma coordinate mesh onto regular coordinates and save as *_rf_itpd.nc')
parser.add_argument('-f', '--fnames-list', nargs='+', default=[])
fnames = parser.parse_args().fnames_list

print('List of files to interpolate:', fnames)

for fname in fnames:
    print('Starting with {}'.format(fname))
    if os.path.exists(fname):

        f = xr.open_dataset(fname)

        for q in ['edir', 'edn', 'eup', 'abso', 'hhl', 'lwc']:
            d = f[q].data[:]
            d = d[:int(d.shape[0]/2),:,:]
            d1 = d[:,:int(d.shape[1]/2),:]
            d2 = d[:,int(d.shape[1]/2):,:][:,::-1,:]
            d = (d1 + d2) / 2
            s = q + '_rf'
            f[s] = ((s+'_dim3',s+'_dim2',s+'_dim1'),d)

        for q in ['edir_srfc', 'edn_srfc', 'eup_srfc', 'h_srfc']:
            d = f[q].data[:]
            d = d[:int(d.shape[0]/2),:]
            d1 = d[:,:int(d.shape[1]/2)]
            d2 = d[:,int(d.shape[1]/2):][:,::-1]
            d = (d1 + d2) / 2
            s = q + '_rf'
            f[s] = ((s+'_dim3',s+'_dim2'),d)
        f.to_netcdf(fname.replace('.nc','_rf.nc'))
        print('Done with {}'.format(fname))
    else:
        print('Could not find {}.'.format(fname))
