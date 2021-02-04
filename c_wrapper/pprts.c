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
#include "f2c_solver_ids.h"

static char help[] = "This is the C wrapper interface to the pprts solver environment.\n\n";

int collapseindex = 1;

int master(int fcomm) {
  int Nx=3, Ny=12, Nz=64;
  int solver_id=SOLVER_ID_PPRTS_3_10;
  double dx=100,dy=100, dz=40.414518843273818;
  float phi0=0, theta0=60;
  float albedo=.2;

  float *hhl   = (float *)malloc((Nz+1) *sizeof(float) );

  float *kabs  = (float *)malloc(Nz*Ny*Nx * sizeof(float) );
  float *ksca  = (float *)malloc(Nz*Ny*Nx * sizeof(float) );
  float *g     = (float *)malloc(Nz*Ny*Nx * sizeof(float) );
  float *planck= (float *)malloc((Nz+1)*Ny*Nx * sizeof(float) );

  float *edir  = (float *)malloc((Nz+1)*Ny*Nx * sizeof(float) );
  float *edn   = (float *)malloc((Nz+1)*Ny*Nx * sizeof(float) );
  float *eup   = (float *)malloc((Nz+1)*Ny*Nx * sizeof(float) );
  float *abso  = (float *)malloc(Nz*Ny*Nx * sizeof(float) );


  for(int j=0;j<Ny;j++) {
    for(int i=0;i<Nx;i++) {
      for(int k=0;k<Nz;k++) {
        int ind = k + Nz*i + Nz*Nx*j; /* index for [Ny][Nx][Nz] */

        kabs   [ind] = 1e-6;
        ksca   [ind] = 1e-6;
        g      [ind] =   .0;

        if(k==0) {
          kabs   [ind] = 1e-8;
        }

        if(k==5 && i==1 && j==1) {
          kabs   [ind] = 1e-30;
          ksca   [ind] = 1e-3;
          g      [ind] =   .0;
        }


      }
    }
  }

  for(int j=0;j<Ny;j++) {
    for(int i=0;i<Nx;i++) {
      for(int k=0;k<Nz+1;k++) {
        int ind = k + (Nz+1)*i + (Nz+1)*Nx*j; /* index for [Ny][Nx][Nz] */

        planck [ind] =  10;/* 5.67e-8*273.*273.*273.*273./3.1415; */
        edir   [ind] =  10;
      }
    }
  }

  hhl[Nz] = 0;
  for(int k=Nz;k>0;k--)
    hhl[k-1] = hhl[k]+dz;

  pprts_f2c_init(fcomm,&solver_id,&Nz,&Nx,&Ny, &dx,&dy, hhl, &phi0, &theta0, &collapseindex);
  pprts_f2c_set_global_optical_properties(Nz,Nx,Ny, &albedo, kabs, ksca, g, planck);
  pprts_f2c_solve(fcomm, 1. );
  pprts_f2c_get_result(Nz,Nx,Ny, edn, eup, abso, edir);

  pprts_f2c_destroy(0);

  for(int j=0;j<Ny/10;j++){
    int i=1;
    printf("\n %d edir ",j);
    for(int k=0;k<Nz+1;k++)
      printf(" %f ",edir[(Nz+1)*Nx*j + (Nz+1)*i + k]);
    printf("\n edn ");
    for(int k=0;k<Nz+1;k++)
      printf(" %f ",edn[(Nz+1)*Nx*j + (Nz+1)*i + k]);
    printf("\n eup ");
    for(int k=0;k<Nz+1;k++)
      printf(" %f ",eup[(Nz+1)*Nx*j + (Nz+1)*i + k]);
    printf("\n abso      ");
    for(int k=0;k<Nz;k++)
      printf(" %f ",abso[Nz*Nx*j + Nz*i + k]*dz);
    printf("\n");
    printf("\n");
  }

  free(hhl   );
  free(kabs  );
  free(ksca  );
  free(g     );
  free(planck);
  free(edir  );
  free(edn   );
  free(eup   );
  free(abso  );

  printf(" finished execution success!\n");
  return 0;
}
int slave(int fcomm) {
  int Nx, Ny, Nz;
  int solver_id=0;
  double dx,dy;
  float phi0, theta0;
  float *albedo = NULL;

  float *kabs   = NULL;
  float *ksca   = NULL;
  float *g      = NULL;
  float *planck = NULL;
  float *hhl    = NULL;

  float *edir   = NULL;
  float *edn    = NULL;
  float *eup    = NULL;
  float *abso   = NULL;

  pprts_f2c_init(fcomm,&solver_id,&Nz,&Nx,&Ny, &dx,&dy, hhl, &phi0, &theta0, &collapseindex);
  pprts_f2c_set_global_optical_properties(Nz,Nx,Ny, albedo, kabs, ksca, g, planck);
  pprts_f2c_solve(fcomm, 1.);
  pprts_f2c_get_result(Nz,Nx,Ny, edn, eup, abso, edir);

  pprts_f2c_destroy(0);
  return 0;
}

int main(int argc, char *argv[]) {
  int        numprocs, myid, fcomm;

  MPI_Init(&argc,&argv);
  PetscInitialize(&argc,&argv,(char*)0,help);
  PetscInitializeFortran();

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
  fprintf(stderr,"I'm number %d of %d.\n", myid, numprocs);

  if(myid==0) {
    master(fcomm);
  }
  else {
    slave(fcomm);
  }

  PetscFinalize();
  MPI_Finalize();
  return(0);

}

