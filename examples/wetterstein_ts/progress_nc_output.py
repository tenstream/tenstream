import numpy as np
import xarray as xr
import os

#for cld in ['cld']:
#    for theta in [25, 40]:
#        for s in ['rayli_ac_sc', 'twostream_sc', 'gomtrc_ac_sc']:
#            fname = s + '_180_{}_{}.nc'.format(theta, cld)
for fname in ['twostream_180_40_cld.nc']:
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
