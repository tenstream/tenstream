import numpy as np
import xarray as xr
from pylab import *
import argparse

parser = argparse.ArgumentParser(description='Plot the results from libRadtran triangular, rectangular and from tenstream.')
parser.add_argument('sza', metavar='sza', type=int, help='the sun zenith angle')
parser.add_argument('phi0', metavar='phi0', type=int, help='phi0')
parser.add_argument('quantity', metavar='quantity', help='The phyical quantity you want to plot.')
args = parser.parse_args()

theta0 = args.sza
phi0 = args.phi0
quantity = args.quantity
print(theta0, phi0, quantity)
elev_file = 'uvspec_elevation'

# libRadtran rectangular results
with xr.open_dataset('libRadtran_rectangular_merged_{}_{}.nc'.format(theta0, phi0)) as data:
  plot(data['y'][:].data, data[quantity].data[:].mean((0,2)), label='rectangular')

# libRadtran triangular results
with xr.open_dataset('libRadtran_triangular_merged_{}_{}.nc'.format(theta0, phi0)) as rad_data:
  with xr.open_dataset('triangular_srfc_file.nc') as srfc_data:
    plot(np.array([ np.array([srfc_data['vertices'][indice] for indice in triangle]).mean(axis=0) for triangle in srfc_data['triangles'] ])[:, 1], rad_data[quantity], '.', label='triangular')

# tenstream results
with xr.open_dataset('out_pprts_hill.nc') as data:
  tenstream_res_edir_srfc = np.average(data['0' + quantity].data[:][:, -1], axis=1)
  plot(np.arange(0, data.Ny * data.dy, data.dy), data['0' + quantity].data[:, :, -1].mean(axis=1), label='tenstream') 

legend()
plot_name = '' + str(theta0) + '_' + str(phi0) + '_' + quantity + '.pdf'
savefig(plot_name)
print('Plot saved to ' + plot_name) 
