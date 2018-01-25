from pylab import *
import netCDF4 as NC
import sys

Nz = int(sys.argv[1])
Ncells = 24

D=NC.Dataset('lwc.nc','w')

D.createDimension('ncells', Ncells)
D.createDimension('hhl_level', Nz+1)
D.createDimension('hhl', Nz)

hhl=D.createVariable('height',float32, dimensions=('hhl',))
for i in range(Nz):
    hhl[i] = 250 + 500*i
hhl[:] = hhl[::-1]

hl=D.createVariable('height_level',float32, dimensions=('hhl_level',))
for i in range(Nz+1):
    hl[i] = 0 + 500*i
hl[:] = hl[::-1]

lwc=D.createVariable('lwc',float32, dimensions=('hhl','ncells'))

lwc[:] = 0
lwc[Nz/2,[18,19,20,21,23]] = 1

D.sync()
D.close()
