/* -------------------------------------------------------------------------*/
/*  This file is part of the tenstream solver.                              */
/*                                                                          */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation, either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License       */
/*  along with this program.  If not, see <http://www.gnu.org/licenses/>.   */
/*                                                                          */
/*  Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>               */
/* -------------------------------------------------------------------------*/

void f2c_tenstream_rrtmg(int fcomm,                   // MPI_Comm_c2f(MPI_COMM_WORLD)
    int *Nz, int *Nx,int *Ny,                         // size of local subdomain
    double *dx, double *dy,                           // horizontal grid spacing in [m]
    double *phi0, double *theta0,                     // Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    double *albedo_thermal, double* albedo_solar,     // broadband ground albedo for solar and thermal spectrum
    char *atm_filename,                               // Filename of background atmosphere file. ASCII file with columns: z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    int *lthermal, int *lsolar,                       // Compute solar or thermal radiative transfer. Or compute both at once.
    int *Nz_merged,                                   // will determine the number of layers of the result
    double **edir, double **edn, double **eup,        // fluxes edir, edn, eup have shape(Nz_merged+1, Nx, Ny)
    double **abso,                                    // abso just (Nz_merged, Nx, Ny)
    double *d_plev,                                   // pressure on layer interfaces    [hPa]
    double *d_tlev,                                   // temperature on layer interfaces [K]
    double *d_lwc,                                    // liq water content               [g/kg]
    double *d_reliq,                                  // effective radius                [micron]
    double *d_iwc,                                    // ice water content               [g/kg]
    double *d_reice,                                  // ice effective radius            [micron]
    int *nprocx,                                      // number of processors in x
    int *nxproc,                                      // local size of subdomain along x, has size(nprocx)
    int *nprocy,                                      // number of processors in y
    int *nyproc                                       // local size of subdomain along y, has size(nprocy)
    );

void f2c_destroy_tenstream_rrtmg(int *lfinalizepetsc);
