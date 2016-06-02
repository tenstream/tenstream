#include <stdio.h>
#include <math.h>
#include "fpda_rrtm_lw.h"

double Planck (double T, double lambda) { 
  const double const_kb=1.38065e-23; // Boltzmann constant [J/K] 
  const double const_h=6.626068e-34; // Planck constant [J*s]
  const double const_c=2.99792458e8; // speed of light [m/s]
  return 2.0*const_h*const_c*const_c/(lambda*lambda*lambda*lambda*lambda)/(exp(const_h*const_c/(lambda*const_kb*T))-1);
}

double Planck_int(double T, double l1, double l2) {
  double integral = 0;
  int Nsamples=(int)ceil((l2-l1)/.1e-6);
  double dlambda=(l2-l1)/Nsamples;
  for(int i=0;i<Nsamples;i++) {
    double lambda = l1 + .5*dlambda + i*dlambda;
    integral += Planck(T, lambda)*dlambda;
   }
  return integral;
}

int main() {
      int nlay=20;
      double plev[nlay+1];
      double tlay[nlay];
      double h2ovmr[nlay]   , o3vmr[nlay]    , co2vmr[nlay]   , ch4vmr[nlay]   , n2ovmr[nlay] , o2vmr[nlay];
      double cfc11vmr[nlay] , cfc12vmr[nlay] , cfc22vmr[nlay] , ccl4vmr[nlay];

      int nbands;
      double * band_lbound;   // [nbands]    
      double * band_ubound;   // [nbands]    
      double * tau_gas;       // [nlay*nbands]       
      double * ck_wght_lw;    // [nlay*nbands]    

      int ib,k;

      double srf_emission,B;

      for(k=0;k<nlay+1;k++) {
        plev[k]=1e3 - (nlay-k)*1e3/(nlay);
        printf("pres(%4i) %g\n",k,plev[k]);
      }
      for(k=0;k<nlay;k++) {
        tlay     [ k ] = 288. - (nlay-k-1)*50/(nlay-1);
        h2ovmr   [ k ] = 0;//9e-6;
        o3vmr    [ k ] = 0;//5e-9;
        co2vmr   [ k ] = 0;//400e-6;
        ch4vmr   [ k ] = 0;//10e-6;
        n2ovmr   [ k ] = 0;//320e-9;
        o2vmr    [ k ] = 0;//.209;
        cfc11vmr [ k ] = 0;
        cfc12vmr [ k ] = 0;
        cfc22vmr [ k ] = 0;
        ccl4vmr  [ k ] = 0;

        printf("lvl(%7i) :: T %7g :: h2o %7g :: o3 %7g :: co2 %7g :: ch4 %7g :: n2o %7g :: o2 %7g \n",k, tlay[k], h2ovmr[k], o3vmr[k], co2vmr[k], ch4vmr[k], n2ovmr[k], o2vmr[k]);
      }

      fpda_rrtm_lw(&nlay,plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
          cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,
          &nbands, &band_lbound,&band_ubound, &ck_wght_lw, &tau_gas );

      k = nlay-1; // srf index is nlay (last layer)

      double tau=0;
      srf_emission=0;
      for( ib=0; ib<nbands; ib++ ) {
        B = Planck_int( tlay[k], 1e-2/band_ubound[ib], 1e-2/band_lbound[ib] );
        tau = tau_gas[ib*nlay + k];
        srf_emission = srf_emission + ck_wght_lw[ib*nlay + k] *B *(1.-exp(-tau));

        printf("Band %4i :: Wavenumber %8g  %8g [cm^-1] :: Wavelength interval %8g  %8g [um] :: ck weight %8e  integrated planck %8g tau %8g\n",k,band_lbound[ib],band_ubound[ib],1e4/band_ubound[ib],1e4/band_lbound[ib], ck_wght_lw[ib*nlay + k], srf_emission,tau);
      }

      printf("\n Nbands %i :: total surface emission %e :: stefan boltzmann %e \n",nbands,srf_emission * 3.1415, 5.67e-8*tlay[k]*tlay[k]*tlay[k]*tlay[k]);


      return 0;
}
