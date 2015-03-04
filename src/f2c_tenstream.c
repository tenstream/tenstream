#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>                                                                                                                    
#include "petscsys.h" 

void tenstr_f2c_init(int fcomm, int *Nx,int *Ny,int *Nz,double *dx,double *dy,float *hhl, float *phi0, float *theta0, float *albedo );
void tenstr_f2c_set_global_optical_properties(int Nx,int Ny,int Nz, float *kabs, float *ksca, float *g, float *planck);
void tenstr_f2c_solve(double edirTOA);
void tenstr_f2c_destroy();
void tenstr_f2c_get_result(int Nx,int Ny,int Nz, float *edir, float *edn, float *eup, float *abso);

static char help[] = "This is the C wrapper interface to the Tenstream solver environment.\n\n";

int master(int fcomm) {
  int    Nx=20, Ny=20, Nz=10;
  double dx=50,dy=50;
  float phi0=180, theta0=0;
  float albedo=.05;

  float *hhl   = (float *)malloc((Nz+1) *sizeof(float) );

  float *kabs  = (float *)malloc(Nz*Ny*Nx * sizeof(float) );
  float *ksca  = (float *)malloc(Nz*Ny*Nx * sizeof(float) );
  float *g     = (float *)malloc(Nz*Ny*Nx * sizeof(float) );
  float *planck= (float *)malloc((Nz+1)*Ny*Nx * sizeof(float) );

  float *edir  = (float *)malloc((Nz+1)*Ny*Nx * sizeof(float) );
  float *edn   = (float *)malloc((Nz+1)*Ny*Nx * sizeof(float) );
  float *eup   = (float *)malloc((Nz+1)*Ny*Nx * sizeof(float) );
  float *abso  = (float *)malloc(Nz*Ny*Nx * sizeof(float) );


  for(int i=0;i<Nx;i++) {
    for(int j=0;j<Ny;j++) {
      for(int k=0;k<Nz;k++) {
        int ind = i + Nx*j + Nx*Ny*k; /* index for [Nz][Ny][Nx] */
        kabs   [ind] = 1e+3;
        ksca   [ind] = 1e-40;
        g      [ind] =   .5;
        planck [ind] =   5.67e-8*273.*273.*273.*273./3.1415;
        edir   [ind] =  -1.;
      }
      int ind = i + Nx*j + Nx*Ny*Nz;
      edir   [ind] =  -1.;
      planck [ind] =  5.67e-8*273.*273.*273.*273./3.1415; /* surface emission */
    }
  }             

  hhl[Nz] = 0;
  for(int k=Nz;k>0;k--) 
    hhl[k-1] = hhl[k]+40.;

  tenstr_f2c_init(fcomm,&Nx,&Ny,&Nz, &dx,&dy, hhl, &phi0, &theta0,&albedo);
  tenstr_f2c_set_global_optical_properties(Nx,Ny,Nz, kabs, ksca, g, planck);
  tenstr_f2c_solve(0.);
  tenstr_f2c_get_result(Nx,Ny,Nz, edir,edn,eup,abso);

  tenstr_f2c_destroy();

  for(int j=0;j<Ny;j++){
    printf("\n edir ");
    for(int k=0;k<Nz+1;k++)
      printf(" %f ",edir[Nx*j + Nx*Ny*k]);
    printf("\n edn ");
    for(int k=0;k<Nz+1;k++)
      printf(" %f ",edn[Nx*j + Nx*Ny*k]);
    printf("\n eup ");
    for(int k=0;k<Nz+1;k++)
      printf(" %f ",eup[Nx*j + Nx*Ny*k]);
    printf("\n abso      ");
    for(int k=0;k<Nz;k++)
      printf(" %f ",abso[Nx*j + Nx*Ny*k]);
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
  int    Nx, Ny, Nz;
  double dx,dy;
  float phi0, theta0;
  float albedo;

  float *kabs   = NULL;
  float *ksca   = NULL;
  float *g      = NULL;
  float *planck = NULL;
  float *hhl    = NULL;

  float *edir   = NULL;
  float *edn    = NULL;
  float *eup    = NULL;
  float *abso   = NULL;

  tenstr_f2c_init(fcomm,&Nx,&Ny,&Nz, &dx,&dy, hhl, &phi0, &theta0,&albedo);
  tenstr_f2c_set_global_optical_properties(Nx,Ny,Nz, kabs, ksca, g, planck);
  tenstr_f2c_solve(1.);
  tenstr_f2c_get_result(Nx,Ny,Nz, edir,edn,eup,abso);

  tenstr_f2c_destroy();
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

