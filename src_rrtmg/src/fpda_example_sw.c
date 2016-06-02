#include <stdio.h>
#include <math.h>
#include "fpda_rrtm_sw.h"

int main() {
      int nlay=20;
      double plev[nlay+1];
      double tlay[nlay];
      double h2ovmr[nlay], o3vmr[nlay], co2vmr[nlay], ch4vmr[nlay], n2ovmr[nlay], o2vmr[nlay];
      double edir  [nlay+1];
      double edir_i[nlay+1];

      int nbands;
      double * band_lbound; // [nbands]    
      double * band_ubound; // [nbands]    
      double * weights;     // [nbands]
      double * tau_gas;     // [nlay*nbands] 
      double * tau_ray;     // [nlay*nbands]

      double tmp=0;

      for(int k=0;k<nlay+1;k++) {
        plev  [k] = k*1e3/(nlay+1);
        edir_i[k] = 0;
      }
      for(int k=0;k<nlay;k++) {
        tlay  [k] = 300. - (nlay-k-1) *50./nlay;
        h2ovmr[k] = 9e-6;
        o3vmr [k] = 5e-9;
        co2vmr[k] = 400e-6;
        ch4vmr[k] = 10e-6;
        n2ovmr[k] = 320e-9;
        o2vmr [k] = .209;
      }

      // Call rrtm routine to get absorption coefficients
      fpda_rrtm_sw(&nlay,plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
          &nbands, &band_lbound,&band_ubound, &weights, &tau_gas, &tau_ray);

      for(int ib=0;ib<nbands;ib++) {
        tmp += weights[ib];
        printf("Band %4i :: Wavenumber %8g  %8g [cm^-1] :: Wavelength interval %8g  %8g [nm] :: ck weight %6.3e  \n",ib,band_lbound[ib],band_ubound[ib],1e7/band_ubound[ib],1e7/band_lbound[ib], weights[ib]);
      }

      // Calculate direct radiation...
      for(int ib=0;ib<nbands;ib++) {
        edir[0] = weights[ib]; // init at TOA with ck weight
        for(int k=0; k<nlay; k++) {
          int ck_ind = ib*nlay + k;                         // index to optical properties in flattened 2d array with fortran ordering
          double dtau = tau_gas[ck_ind] + tau_ray[ck_ind];  // direct radiation is subject to extinction due to gas absorption and rayleigh scattering
          edir[k+1] = edir[k] * exp(-dtau);
        }

        for(int k=0; k<nlay+1; k++)
          edir_i[k] += edir[k];         // add subband edir to spectral integral

      }

      for(int k=0; k<nlay+1; k++) 
        printf("direct radiation(%4i) :: %g \n",k,edir_i[k]);

      printf("\n Bands %i :: total solar weights %g \n\n",nbands,tmp);

      return 0;
}
