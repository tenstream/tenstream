import numpy as np
import xarray as xr
from scipy.interpolate import griddata
import os
import glob

### Generating uvspec input files for calculations
def run_uvspec_job(datadict):
    """
    Generate a LibRadtran input with a template and run it.
    Input file is a Temporary file and
    Output is saved in a temporary folder
    """
    import os
    import shutil
    import subprocess
    import tempfile
    from jinja2 import Template

    inp_template = Template(
    """
data_files_path {{data_files_path}}
atmosphere_file {{atmosphere_file}}
source solar    {{source_solar}} per_nm

{% if output_process %}
output_process  {{output_process}}
{% else %}
output_process  sum
{% endif %}

mol_abs_param   {{mol_abs_param}}
wavelength_index {{wavelength_index}}

albedo          {{albedo}}
sza             {{sza}}
phi0            {{phi0}}

{% if aerosol %}
{{aerosol}}
{% endif %}

{% if mc_forward_output %}
mc_forward_output {{mc_forward_output}}
{% endif %}

{% if mc_backward_output %}
mc_backward
mc_backward_output {{mc_backward_output}}
{% endif %}

{% if zout %}
zout            {{zout}}
{% endif %}

rte_solver      {{rte_solver}}
{% if ipa_3d %}
ipa_3d
{% endif %}

{% if mc_tenstream %}
mc_tenstream    {{mc_tenstream}}
{% endif %}

{% if mc_ipa %}
mc_ipa
{% endif %}

{% if mc_vroom %}
mc_vroom {{mc_vroom}}
{% endif %}

{% if mc_delta_scaling %}
mc_delta_scaling {{mc_delta_scaling}}
{% endif %}

{% if mc_photons %}
mc_photons      {{mc_photons}}
mc_minphotons   10
{% endif %}

mc_sample_grid  {{mc_sample_grid}}

{% if wc_file_3D %}
wc_file 3D      {{wc_file_3D}}
wc_properties   hu interpolate
{% endif %}

{% if ic_file_3D %}
ic_file 3D      {{ic_file_3D}}
ic_properties   yang interpolate
{% endif %}

{% if outdir %}
mc_basename     {{outdir}}/mc
{% endif %}

{% if mc_elevation_file %}
mc_elevation_file {{mc_elevation_file}}
{% endif %}

{% if mc_surfaceparallel %}
mc_surfaceparallel
{% endif %}

#quiet
verbose
    """)
    if 'outdir' in datadict:
        with open('uvspec_'+datadict['outdir']+'.inp', 'w') as fh:
            fh.write(inp_template.render(**datadict))
            fh.flush()

            if not os.path.exists(datadict['outdir']):
                os.mkdir(datadict['outdir'])
            else:
                if os.path.exists( os.path.join(datadict['outdir'], 'mc.flx.spc' ) ):
                    print("Skipping job for {} because output already exists".format(datadict['outdir']))
                    return
                else:
                    print("Resuming job for {} ".format(datadict['outdir']))

            ts_opt = ''
            if 'tenstream_options' in datadict:
                ts_opt = ' '.join(map(str, datadict.pop('tenstream_options')))

            if 'mc_tenstream' in datadict:
                cmd = ['salloc', '-n', '64', '--time=08:00:00', '--mem-per-cpu=4G',  '-C', 'GPU', '-p', 'vis,cluster,ws', 'bash', '-c', 'mpirun uvspec_mpi {} {}'.format(fh.name, ts_opt)]
            else:
                cmd = ['srun', '-n', '1', '--time=08:00:00', '--mem=1G']
                cmd += [os.path.join(libRadtran, 'bin', 'uvspec_mpi'), fh.name]

            print("calling subprocess", cmd)
            sp = subprocess.call(cmd)
            [ load_flx_spc(flxfile) for flxfile in glob.glob(os.path.join(datadict['outdir'],'*.flx.spc')) ]
            [ load_abs_spc(absfile) for absfile in glob.glob(os.path.join(datadict['outdir'],'*.abs.spc')) ]
    else:
        print('outdir not in datadict, not sure where to put the results, skipping...')

def wrf_HGT_2_elevation(hgt, dx, dy, fill_value='nearest', max_resolution_decimals=1):
    """
    hgt has to be a 2D array [km]
    dx, dy are scalars of grid resolution [km]

    fill_value:
      'nearest' : sets extrapolated heigths to nearest neighbor height
      'lowest' : sets extrapolated heigths to globally lowest height value
    """
    Ny, Nx = hgt.shape
    y, x = np.meshgrid(np.arange(Ny)*dy, np.arange(Nx)*dx, indexing='ij')

    y_stag, x_stag = np.meshgrid(np.arange(-.5, Ny)*dy, np.arange(-.5, Nx)*dx, indexing='ij')

    def append_border_rows_and_cols(A):
        """ append rows and cols around the border, so that we get periodic boundary conditions """
        A = np.vstack((A[:,0], A.T, A[:,-1])).T
        A = np.vstack((A[0,:].T, A, A[-1,:].T))
        return A

    periodic_hgt = append_border_rows_and_cols(hgt)
    y, x = np.meshgrid(np.arange(-1,Ny+1)*dy, np.arange(-1,Nx+1)*dx, indexing='ij')
    new_hgt = griddata((y.flatten(),x.flatten()), periodic_hgt.flatten(), (y_stag.flatten(), x_stag.flatten()), method='linear').reshape(np.shape(y_stag))

    if fill_value == 'lowest':
        new_hgt[np.isfinite(new_hgt)==False] = np.min(hgt)
    elif fill_value == 'nearest':
        nearest_hgt = griddata((y.flatten(),x.flatten()), periodic_hgt.flatten(), (y_stag.flatten(), x_stag.flatten()), method='nearest').reshape(np.shape(y_stag))
        new_hgt[np.isfinite(new_hgt)==False] = nearest_hgt[np.isfinite(new_hgt)==False]
    else:
        raise('Unsupported fill_value method')

    new_hgt = np.round(new_hgt, decimals= max_resolution_decimals)

    return new_hgt


def tenstr_hill_2_elevation_file(fname, outfname=None, xslice=None, yslice=None):
    with xr.open_dataset(fname) as D:
        dx, dy = D.dx, D.dy
        hgt = (D['0zt'].max() - D['0zt'][:,:,-1]).data
        if xslice is None:
            xslice = slice(None)
        if yslice is None:
            yslice = slice(None)
        hgt = hgt[yslice,xslice]
        new_hgt = wrf_HGT_2_elevation(hgt, dx, dy, max_resolution_decimals=3)

    if outfname is None:
        outfname='elevation.{}'.format(os.path.basename(fname))

    D = xr.Dataset({
        'elevation': xr.DataArray(new_hgt*1e-3, dims=('ny','nx')),
        'dx': xr.DataArray(dx*1e-3),
        'dy': xr.DataArray(dy*1e-3),
        })

    D.to_netcdf(outfname)

def tenstr_hill_2_cloud_file(fname, outfname=None, max_vert_resolution_decimals=-1, xslice=None, yslice=None):
    with xr.open_dataset(fname) as D:
        if xslice is None:
            xslice = slice(None)
        if yslice is None:
            yslice = slice(None)
        qcloud = D['0lwc'][yslice, xslice, :].data.T * 1e3
        qcloud[:,0,0][qcloud[:,0,0]<1e-10] = 1e-10
        zslice  = slice(None, qcloud.shape[0])
        zslice1 = slice(None, qcloud.shape[0]+1)

        hhl = D['0zt'][yslice, xslice, zslice1].data.T

        dx, dy = D.dx, D.dy

    new_hhl = np.unique(np.round(hhl, decimals=max_vert_resolution_decimals))

    cld_hgt = (hhl[1:] + hhl[1:])/2
    new_cld_hgt = (new_hhl[1:] + new_hhl[1:])/2

    Nz, Nx, Ny = qcloud.shape

    _, x, y = np.meshgrid(np.arange(Nz), np.arange(Nx)*dx, np.arange(Ny)*dy, indexing='ij')
    i_pts = tuple([ _.flatten() for _ in (cld_hgt,x,y) ])


    newz, newx, newy = np.meshgrid(new_cld_hgt, np.arange(Nx)*dx, np.arange(Ny)*dy, indexing='ij')
    o_pts = tuple([ _.flatten() for _ in (newz,newx,newy) ])

    new_qcloud = griddata(i_pts, qcloud.flatten(), o_pts, method='nearest').reshape(newz.shape)

    reff = np.zeros_like(new_qcloud)
    reff[new_qcloud>0] = 10

    if outfname is None:
        outfname='wcloud.{}.{}'.format(it, os.path.basename(fname))

    print("max LWC {} [g/kg] shape({})".format(np.max(new_qcloud), new_qcloud.shape))

    D = xr.Dataset({
        'lwc': xr.DataArray(new_qcloud.T, dims=('ny','nx','nz')),
        'reff': xr.DataArray(reff.T, dims=('ny','nx','nz')),
        'z': xr.DataArray(new_hhl*1e-3, dims=('nz_lev')), },
        attrs={
            'cldproperties': 3,
            'dx': dx*1e-3,
            'dy': dy*1e-3,
        })
    D.to_netcdf(outfname)


def load_flx_spc(fname,d=dict(),ret=True,delete=False):
        from numpy import loadtxt,arange,zeros,array,save,load
        import os

        if os.path.splitext(fname) == '.npy':
                if not ret:
                	return
                returnarr = load(fname+'.npy')
        else:
                print('\t \t ... converting flux file',fname)
                x,y,z,edir,edn,eup,uavgdir,uavgdn ,uavgup = loadtxt(fname,unpack=True,usecols=(1,2,3,4,5,6,7,8,9) )
                edir    = edir.reshape    (  ( -1 , int(max(x))+1 , int(max(y))+1 ))
                edn     = edn.reshape     (  ( -1 , int(max(x))+1 , int(max(y))+1 ))
                eup     = eup.reshape     (  ( -1 , int(max(x))+1 , int(max(y))+1 ))
                uavgdir = uavgdir.reshape (  ( -1 , int(max(x))+1 , int(max(y))+1 ))
                uavgdn  = uavgdn.reshape  (  ( -1 , int(max(x))+1 , int(max(y))+1 ))
                uavgup  = uavgup.reshape  (  ( -1 , int(max(x))+1 , int(max(y))+1 ))
                returnarr = array([edir,edn,eup ,uavgdir,uavgdn,uavgup]).swapaxes(2,3)
                save(fname,returnarr)
                if os.path.exists(fname+'.npy'):
                    if delete:
                        os.remove(fname)

        d[str(len(d))+':'+(fname.rsplit('/')[-1]).replace('.out.flx.spc','')] = returnarr
        return d


def load_abs_spc(fname,d=dict(),ret=True, delete=False):
    from numpy import loadtxt,arange,zeros,array,save,load
    import os
    if os.path.splitext(fname) == '.npy':
        if not ret:
            return
        returnarr = load(fname+'.npy')
    else:
        print('\t \t ... converting flux file',fname)
        wvl,x,y,z,abs = loadtxt(fname,unpack=True )
        abs   = abs.reshape((-1 , int(max(x))+1 , int(max(y))+1 ))
        returnarr = abs.swapaxes(1,2)
        save(fname,returnarr)
        if os.path.exists(fname+'.npy'):
            if delete:
                os.remove(fname)

    d[str(len(d))+':'+(fname.rsplit('/')[-1]).replace('.out.flx.spc','')] = returnarr
    return d

# HILL JOB

import os
import multiprocessing

homedir = os.environ['HOME']
if 'LIBRADTRAN_PATH' in os.environ:
    libRadtran = os.environ['LIBRADTRAN_PATH']
elif os.path.exists(os.path.join(homedir, 'libRadtran')):
    libRadtran = os.path.join(homedir, 'libRadtran')
else:
    raise RuntimeError('Could not find libRadtran')

tenstr_hill_fname='out_pprts_hill.nc'
elev_file = 'uvspec_elevation.nc'
wcloud_file = 'wcloud_tenstr_hill.nc'

tenstr_hill_2_elevation_file(tenstr_hill_fname, elev_file)
tenstr_hill_2_cloud_file(tenstr_hill_fname, wcloud_file)

with xr.open_dataset(tenstr_hill_fname) as D:
    outdir_srfc = 'out_hill_{}_{}_srfc'.format(int(D.theta0), int(D.phi0))
    outdir_atm  = 'out_hill_{}_{}_atm'.format(int(D.theta0), int(D.phi0))
    datadict = {
	    'data_files_path': os.path.join(libRadtran, 'data'),
	    'atmosphere_file': os.path.join(libRadtran, 'data', 'atmmod', 'afglus.dat'),
	    'source_solar': os.path.join(libRadtran, 'data', 'atlas_plus_modtran'),
	    'mol_abs_param': 'kato2',
	    'wavelength_index': '1 32',
	    'albedo': D.Ag_solar,
	    'sza': D.theta0,
	    'phi0': D.phi0+180, # 0 is South

	    'rte_solver': 'montecarlo',

	    'mc_sample_grid': '{} {} {} {}'.format(D.Nx, D.Ny, D.dx*1e-3, D.dy*1e-3),
	    'wc_file_3D': wcloud_file,
	    'mc_forward_output': 'absorption K_per_day',

	    'mc_photons': '{}'.format(int(1e7)),
	    'outdir': outdir_srfc,

            'mc_elevation_file': elev_file,
            'mc_surfaceparallel': True,
	    }

    joblist = []
    for it in range(20):
        d = datadict.copy()
        d['outdir'] = outdir_srfc+'.iter{}'.format(it)
        joblist.append(d)


    with xr.open_dataset(wcloud_file) as Dcld:
        datadict['zout']   = ' '.join(map(str,Dcld.z.data))

    for it in range(20):
        d = datadict.copy()
        d['outdir'] = outdir_atm+'.iter{}'.format(it)
        joblist.append(d)

    print(joblist)

    with multiprocessing.Pool(40) as p:
        p.map(run_uvspec_job, joblist)

fl_iter = [ np.load(_) for _ in glob.glob(os.path.join(outdir_srfc+'.iter*','*.flx.spc.npy')) ]
fl = np.average(fl_iter, axis=0)
from pylab import *
plot(np.average(fl[0,0,:,:], axis=-1))
axvline((fl.shape[-2]+1)//2, linestyle='--', color='.5')
waitforbuttonpress()
