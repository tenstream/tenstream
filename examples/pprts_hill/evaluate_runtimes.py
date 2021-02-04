import numpy as np
import glob
import os

times_dict = {}
for string in ['_elevation_', '_triangular_']:
    times = np.array([ np.loadtxt(_) for _ in glob.glob(os.path.join('out_hill_0_0_srfc' + string + '.iter*', 'execution_time.txt')) ]).flatten
    times_dict[string] = np.amin(times, axis=0)
np.savetxt('t_tri_over_t_el.txt', np.array([dict['_triangular_'] / dict['_elevation_']]), fmt='%.2f')
