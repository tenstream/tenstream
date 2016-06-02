#ifndef _rrtm_sw_h
#define _rrtm_sw_h 1

#if defined (__cplusplus)
extern "C" {
#endif

#include <stdlib.h>

void tenstr_rrtm_sw(
    // INPUT: profiles start at upper boundary
    int *nlay,               // Number of atmospheric layers
    double *plev,            // pressure on interface                       [ nlay+1 ] [hPa]
    double *tlay,            // temperature in layers                       [ nlay   ] [K  ]
    double *h2ovmr,          // H2O volume mixing ratio                     [ nlay   ]
    double *o3vmr,           // O3 volume mixing ratio                      [ nlay   ]
    double *co2vmr,          // CO2 volume mixing ratio                     [ nlay   ]
    double *ch4vmr,          // Methane volume mixing ratio                 [ nlay   ]
    double *n2ovmr,          // Nitrous oxide volume mixing ratio           [ nlay   ]
    double *o2vmr,           // Oxygen volume mixing ratio                  [ nlay   ]

    // OUTPUT:
    // Do not free those arrays in "C" !  Those pointers merely point to an Fortran allocatable
    //
    int    *nbands,          // number of spectral bands including subbands
    double **band_lbound,    // lower bound of wavenumber                   [ nbands ]
    double **band_ubound,    // upper bound of wavenumber                   [ nbands ]
    double **weights,        // correlated-k weigth                         [ nbands ]
    double **tau_gas,        // gas absorption opt depth                    [ nlay,nbands ] ! Fortran ordered ! i.e. first nlay entries are for first band
    double **tau_ray);       // rayleigh opt depth                          [ nlay,nbands ] ! Fortran ordered ! i.e. first nlay entries are for first band




// Same wrapper as above but with 2D arrays in C-style
void ctenstr_rrtm_sw(
    // INPUT: profiles start at upper boundary
    int nlay,                // Number of atmospheric layers
    double *plev,            // pressure on interface                       [ nlay+1 ] [hPa]
    double *tlay,            // temperature in layers                       [ nlay   ] [K  ]
    double *h2ovmr,          // H2O volume mixing ratio                     [ nlay   ]
    double *o3vmr,           // O3 volume mixing ratio                      [ nlay   ]
    double *co2vmr,          // CO2 volume mixing ratio                     [ nlay   ]
    double *ch4vmr,          // Methane volume mixing ratio                 [ nlay   ]
    double *n2ovmr,          // Nitrous oxide volume mixing ratio           [ nlay   ]
    double *o2vmr,           // Oxygen volume mixing ratio                  [ nlay   ]

    // OUTPUT:
    // Do not free those arrays in "C" !  Those pointers merely point to an Fortran allocatable
    //
    int    *nbands,          // number of spectral bands including subbands
    double **band_lbound,    // lower bound of wavenumber           [ nbands ]
    double **band_ubound,    // upper bound of wavenumber           [ nbands ]
    double **wgt,            // planck fractions, total norm is 14  [ nbands ]
    double ***dtaumol,       // gas absorption opt depth            [ nbands ] [nlay] 
    double ***dtauray)       // Rayleigh scattering opt depth       [ nbands ] [nlay] 
{
  static int first=1;
  static double *tau_gas=NULL, *tau_ray=NULL;
  int iv=0;
  
  tenstr_rrtm_sw (&nlay,plev,tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
		nbands, band_lbound, band_ubound, wgt, &tau_gas, &tau_ray);

  if (first) {   /* allocate memory for 2D result fields */
    *dtaumol = calloc (*nbands, sizeof(double *));
    *dtauray = calloc (*nbands, sizeof(double *));
    first = 0;
  }

  for (iv=0; iv<*nbands; iv++) {
    (*dtaumol)[iv] = tau_gas + iv * nlay;
    (*dtauray)[iv] = tau_ray + iv * nlay;
  }
}

  
#if defined (__cplusplus)
}
#endif

#endif

