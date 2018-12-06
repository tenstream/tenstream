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

#include <stdio.h>
#include <stdlib.h>
#include <petscsys.h>
#include <mpi.h>

#include <f2c_pprts.h>

static char help[] = "This is the C wrapper interface to the pprts solver calling the RRTMG routines.\n\n";

int main(int argc, char *argv[]) {
  int        numprocs, myid, fcomm;

  int    Nx=3, Ny=3, Nz=2;
  double dx=500,dy=500, dz=250;
  double phi0=0, theta0=60;
  double albedo_th=.005, albedo_sol=0.2;
  char   atm_filename[] = "afglus.dat";
  int    lsolar=1, lthermal=1;

  int    Nz_merged;
  double *edir, *edn, *eup, *abso;  // will have shape(Nz_merged*Nx*Ny)

  double *d_tlev = malloc(Nx*Ny*(Nz+1) * sizeof(double));
  double *d_plev = malloc(Nx*Ny*(Nz+1) * sizeof(double));

  double *d_lwc   = malloc(Nx*Ny*Nz * sizeof(double));
  double *d_reliq = malloc(Nx*Ny*Nz * sizeof(double));
  double *d_iwc   = malloc(Nx*Ny*Nz * sizeof(double));
  double *d_reice = malloc(Nx*Ny*Nz * sizeof(double));

  MPI_Init(&argc,&argv);
  PetscInitialize(&argc,&argv,(char*)0,help);
  PetscInitializeFortran();

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);

  // Poor mans domain decompostion, we just use 1 row of processes along the x axis and numprocs along y axis
  int nprocx = 1;
  int *nxproc = malloc(sizeof(int));
  nxproc[0] = Nx;

  int nprocy = numprocs;
  int *nyproc = malloc(nprocy*sizeof(int));
  for(int j=0; j<nprocy; j++) {
    nyproc[j] = Ny;
  }


  // Initialize temperature and pressure profile
  for (int k=0; k<Nz+1; k++) {
    for (int i=0; i<Nx; i++) {
      for (int j=0; j<Ny; j++) {
        int ind = j*Nx*(Nz+1) + i*(Nz+1) + k;     // `ind` gives index [k,i,j] (Fortran Order)
        d_plev[ind] = 1013 - k * 200/(Nz+1);  // go from srfc pressure, 200 hPa down
        d_tlev[ind] = 288 - k * 10/(Nz+1);    // and reduce temperature by 10 K
      }
    }
  }

  // and juice it up with a thick cloud
  for (int k=0; k<Nz; k++) {
    for (int i=0; i<Nx; i++) {
      for (int j=0; j<Ny; j++) {
        int ind = j*Nx*Nz + i*Nz + k;     // `ind` gives index [k,i,j] (Fortran Order)
        d_lwc  [ind] = .1;
        d_reliq[ind] = 10;
        d_iwc  [ind] =  0;
        d_reice[ind] = 10;
      }
    }
  }

  // call rrtmg and pprts
  f2c_pprts_rrtmg(fcomm, &Nz, &Nx, &Ny, &dx, &dy, &phi0, &theta0,
      &albedo_th, &albedo_sol, atm_filename, &lthermal, &lsolar,
      &Nz_merged, &edir, &edn, &eup, &abso, d_plev, d_tlev,
      d_lwc, d_reliq, d_iwc, d_reice, &nprocx, nxproc, &nprocy, nyproc);

  // and the result will be reachable again in Fortran Order
  if(edir) {
      fprintf(stdout, " rank :: level :: edir edn eup abso\n");
  } else {
      fprintf(stdout, " rank :: level :: edn eup abso\n");
  }
  //for (int i=0; i<Nx; i++) {
  //for (int j=0; j<Ny; j++) {
  for (int k=0; k<Nz_merged; k++) {
      //int ind = k + i*Nz_merged + j*Nx*Nz_merged;  // `ind` gives index [k,i,j] (Fortran Order)
      //fprintf(stdout, "%d :: %d %d %d :: %d :: %f \n", myid, k, i, j, ind, edir[ind]);
      if(edir) {
          fprintf(stdout, "%d :: %d :: %f %f %f %f\n", myid, k, edir[k], edn[k], eup[k], abso[k]);
      } else {
          fprintf(stdout, "%d :: %d :: %f %f %f\n", myid, k, edn[k], eup[k], abso[k]);
      }
  }
  int k = Nz_merged;
  if(edir) {
      fprintf(stdout, "%d :: %d :: %f %f %f %f\n", myid, k, edir[k], edn[k], eup[k], abso[k]);
  } else {
      fprintf(stdout, "%d :: %d :: %f %f %f\n", myid, k, edn[k], eup[k], abso[k]);
  }
  //}
  //}

  int lfinalizepetsc = 0; // dont drop the Petsc environment if we called initialize here in the C program
  f2c_destroy_pprts_rrtmg(&lfinalizepetsc); // deletes the state of the pprts solver

  PetscFinalize(); // have to do it on our own...
  MPI_Finalize();
  return(0);

}

