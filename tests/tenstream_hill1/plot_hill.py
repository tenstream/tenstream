from pylab import *
import netCDF4 as NC


def create_linear_hill():
    from scipy.ndimage.filters import gaussian_filter as gf
    from scipy.ndimage.filters import uniform_filter as uf

    Nx,Ny = 32,64

    A=np.zeros((Nx,Ny))
    A[Nx/2-1,Ny/2-1]=1
    g = uf(A,5)
    g = uf(g,5)
    g = g/np.max(g)

    z,p,_,_,_,_,_,_,_ = np.loadtxt('afglus_100m.dat').T
    Nz = np.shape(z)[-1]

    from scipy.interpolate import interp1d
    fp = interp1d(z,p)

    hhl = np.tile(z,Nx*Ny).reshape((Nx,Ny,Nz))[:,:,::-1]  # hhl now begins at surface 

    hill = hhl.copy()

    for k in range(Nz):
        hill[:,:,k] = hhl[:,:,k] + .5 * g * exp(-hhl[:,:,k]/2)

    for i in range(Nx):
        hill[i,:,:] = hill[Nx/2-1,:,:]

    p = fp(hill)[:,:,::-1]

    D=NC.Dataset('hill_input.nc','w')

    D.createDimension('Nx', Nx)
    D.createDimension('Ny', Ny)
    D.createDimension('Nz', Nz)

    x = D.createVariable('Nx', 'd', ('Nx',))
    y = D.createVariable('Ny', 'd', ('Ny',))
    x[:] = linspace(0,Nx*.1,Nx)
    y[:] = linspace(0,Ny*.1,Ny)

    data = D.createVariable('hhl', 'd', ('Ny', 'Nx', 'Nz'))
    data[:] = np.swapaxes(hill,0,1)

    data = D.createVariable('plev', 'd', ('Ny', 'Nx', 'Nz'))
    data[:] = np.swapaxes(p,0,1)

    D.close()


def create_gauss_hill():
    from scipy.ndimage.filters import gaussian_filter as gf

    Nx,Ny = 32,64

    A=np.zeros((Nx,Ny))
    A[Nx/2-1,Ny/2-1]=1
    g = gf(A,10)
    g = g/np.max(g)

    z,p,_,_,_,_,_,_,_ = np.loadtxt('afglus_100m.dat').T
    Nz = np.shape(z)[-1]

    from scipy.interpolate import interp1d
    fp = interp1d(z,p)

    hhl = np.tile(z,Nx*Ny).reshape((Nx,Ny,Nz))[:,:,::-1]  # hhl now begins at surface 
    hill = hhl.copy()

    for k in range(Nz):
        hill[:,:,k] = hhl[:,:,k] + .5 * g * exp(-hhl[:,:,k]/2)

    for i in range(Nx):
        hill[i,:,:] = hill[Nx/2-1,:,:]

    p = fp(hill)[:,:,::-1]

    D=NC.Dataset('hill_input.nc','w')

    D.createDimension('Nx', Nx)
    D.createDimension('Ny', Ny)
    D.createDimension('Nz', Nz)

    x = D.createVariable('Nx', 'd', ('Nx',))
    y = D.createVariable('Ny', 'd', ('Ny',))
    x[:] = linspace(0,Nx*.1,Nx)
    y[:] = linspace(0,Ny*.1,Ny)

    data = D.createVariable('hhl', 'd', ('Ny', 'Nx', 'Nz'))
    data[:] = np.swapaxes(hill,0,1)

    data = D.createVariable('plev', 'd', ('Ny', 'Nx', 'Nz'))
    data[:] = np.swapaxes(p,0,1)

    D.close()

def create_srtm_hill():
    from scipy.ndimage.filters import gaussian_filter as gf
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open('/home/users/jakub/workspace/srtm/srtm_39_03.tif')
    elev = ds.ReadAsArray()
    wetterstein = (elev[3000:3200, 1100:1450]/1e3).T

    
    Nx,Ny = shape(wetterstein)

    afglus = np.loadtxt('afglus_100m.dat')
    z = afglus[:,0]
    Nz = np.shape(z)[-1]

    p = afglus[:,1]

    from scipy.interpolate import interp1d
    fp = interp1d(z,p,kind='linear', fill_value='extrapolate')

    hhl = np.tile(z,Nx*Ny).reshape((Nx,Ny,Nz))[:,:,::-1]  # hhl now begins at surface 
    hill = hhl.copy()

    for k in range(Nz):
        hill[:,:,k] = hhl[:,:,k] + wetterstein * exp(-hhl[:,:,k]/4.)

    #for i in range(Ny):
    #    hill[:,i,:] = hill[:,Ny/2-1,:]

    def create_var(varname, axis, vals, newaxis, dim_prefix):
        import os
        fp = interp1d(axis, vals, kind='linear', fill_value='extrapolate')
        var = fp(newaxis)[:,:,::-1]

        fname = 'hill_input.nc'
        if os.path.exists(fname):
            D=NC.Dataset(fname,'a')
        else:
            D=NC.Dataset(fname,'w')


        Nx,Ny,Nz = np.shape(var)

        try:
            D.createDimension(dim_prefix+'Nx', Nx)
            D.createDimension(dim_prefix+'Ny', Ny)
            D.createDimension(dim_prefix+'Nz', Nz)

            x = D.createVariable(dim_prefix+'Nx', 'd', (dim_prefix+'Nx',))
            y = D.createVariable(dim_prefix+'Ny', 'd', (dim_prefix+'Ny',))
            x[:] = linspace(0,Nx*.1,Nx)
            y[:] = linspace(0,Ny*.1,Ny)
        except:
            pass

        data = D.createVariable(varname, 'd', (dim_prefix+'Ny', dim_prefix+'Nx', dim_prefix+'Nz'))
        data[:] = np.swapaxes(var,0,1)

        D.close()


    create_var('plev'  , afglus[:,0], afglus[:,1], hill,'lev')

    create_var('tlay',   afglus[:,0], afglus[:,2], (hill[:,:,1:]+hill[:,:,:-1])*.5,'lay')
    create_var('air',    afglus[:,0], afglus[:,3], (hill[:,:,1:]+hill[:,:,:-1])*.5,'lay')
    create_var('o3vmr',  afglus[:,0], afglus[:,4], (hill[:,:,1:]+hill[:,:,:-1])*.5,'lay')
    create_var('o2vmr',  afglus[:,0], afglus[:,5], (hill[:,:,1:]+hill[:,:,:-1])*.5,'lay')
    create_var('h2ovmr', afglus[:,0], afglus[:,6], (hill[:,:,1:]+hill[:,:,:-1])*.5,'lay')
    create_var('co2vmr', afglus[:,0], afglus[:,7], (hill[:,:,1:]+hill[:,:,:-1])*.5,'lay')
    create_var('n2ovmr', afglus[:,0], afglus[:,8], (hill[:,:,1:]+hill[:,:,:-1])*.5,'lay')


#    p = fp(hill)[:,:,::-1]
#
#    D=NC.Dataset('hill_input.nc','w')
#
#    try:
#        D.createDimension('Nx', Nx)
#        D.createDimension('Ny', Ny)
#        D.createDimension('Nz', Nz)
#
#        x = D.createVariable('Nx', 'd', ('Nx',))
#        y = D.createVariable('Ny', 'd', ('Ny',))
#        x[:] = linspace(0,Nx*.1,Nx)
#        y[:] = linspace(0,Ny*.1,Ny)
#    except:
#        pass
#
#    data = D.createVariable('hhl', 'd', ('Ny', 'Nx', 'Nz'))
#    data[:] = np.swapaxes(hill,0,1)
#
#    data = D.createVariable('plev', 'd', ('Ny', 'Nx', 'Nz'))
#    data[:] = np.swapaxes(p,0,1)
#
#    D.close()

#create_gauss_hill()
create_srtm_hill()

D=NC.Dataset('output.nc')

hhl = D['hhl'][:]
#hhl = (hhl.T + np.max(hhl) - hhl[:,:,0]).T
h = np.rollaxis(hhl,2)
h = h + np.max(h) - h[0]

edir = D['edir'][:]
ix = np.shape(edir)[-1]/2-1

print ix, shape(edir), shape(h)

X=np.tile(np.arange(np.shape(edir)[-2])*400,len(edir)).reshape((len(edir),-1))
Y=h[::-1][:len(edir),:,ix]

subplot(211)

contourf(X, Y,edir[:,:,ix])
colorbar()

subplot(212)
[ plot( edir[0,:,ix], 'o-', label=str(ix)) for ix in arange(0,8) ]
grid()
legend()

savefig('hill.pdf')

