import xarray as xr
import numpy as np

Nx = 32
Ny = 16
Nz = 10

WALL=1
ROOF=2

hy = Ny//2 -1
hx = Nx//2 -1
q = Nx//3 -1

domain = np.zeros((Nx, Ny, Nz))

domain[q-1:q+2, hy-1:hy+2, :5] = WALL
domain[q-1:q+2, hy-1:hy+2, 5] = ROOF
domain[q, hy, 2:] = 0
domain[q, hy, 1] = ROOF

domain[hx+q-1:hx+q+2, hy-1:hy+2, :5] = WALL
domain[hx+q-1:hx+q+2, hy-1:hy+2, 5] = ROOF
domain[hx+q, hy, 5] = WALL
domain[hx+q, hy, 5+1] = ROOF

domain[q+2:hx+q-1, hy, 3] = WALL
domain[q+2:hx+q-1, hy, 4] = ROOF



albedo = np.zeros((Nx, Ny, Nz, 6)) - 1
temperature = np.zeros((Nx, Ny, Nz, 6)) - 1

for i,j,k in zip(*np.where(domain == WALL)):
    albedo[i,j,k,:] = .4
    temperature[i,j,k,:] = 287.15
for i,j,k in zip(*np.where(domain == ROOF)):
    albedo[i,j,k,:] = .2
    temperature[i,j,k,:] = 287.15

D = xr.Dataset(
        data_vars=dict(
            albedo = (["x", "y", "z", "face"], albedo),
            temperature = (["x", "y", "z", "face"], temperature),
            ),
        coords=dict(),
        attrs=dict(description="test dataset to compute tenstream RT with buildings"),
        )
D.to_netcdf('two_towers_with_bridge.nc')
