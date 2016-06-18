from pylab import *
import netCDF4 as NC

def create_hill():
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

create_hill()

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

