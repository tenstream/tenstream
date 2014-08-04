#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>                                                                                                                    
#include "petscsys.h" 

void f2c_init_tenstream(int fcomm, int Nx,int Ny,int Nz,double dx,double dy,double *hhl, double phi0, double theta0, double albedo );
void f2c_set_optical_properties(int Nx,int Ny,int Nz, double ***kabs, double ***ksca, double ***g);
void f2c_solve_tenstream(double edirTOA);
void f2c_destroy_tenstream();
void get_result(int Nx,int Ny,int Nz, double ***edir, double ***edn, double ***eup, double ***abso);

static char help[] = "This is the C wrapper interface to the Tenstream solver environment.\n\n";

int main(int argc, char *argv[]) { 

  int Nx=1, Ny=1, Nz=10;
  double dx=70,dy=70;
  double phi0=180, theta0=0;
  double albedo=.05;

  double kabs[Nz][Ny][Nx];
  double ksca[Nz][Ny][Nx];
  double g   [Nz][Ny][Nx];
  double hhl[Nz+1];

  double edir[Nz+1][Ny][Nx];
  double edn [Nz+1][Ny][Nx];
  double eup [Nz+1][Ny][Nx];
  double abso[Nz  ][Ny][Nx];

  int        numprocs, myid, fcomm;

  MPI_Init(&argc,&argv);
  PetscInitialize(&argc,&argv,(char*)0,help);
  PetscInitializeFortran();

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
  fprintf(stderr,"I'm number %d of %d.\n", myid, numprocs);

  if(myid==0) {

    for(int i=0;i<Nx;i++) {
      for(int j=0;j<Ny;j++) {
        for(int k=0;k<Nz;k++) {
          kabs[k][j][i] = 1e-3;
          ksca[k][j][i] = 1e-40;
          g   [k][j][i] =   .5;
          edir[k][j][i] =  -1.;
        }
        edir[Nz][j][i] =  -1.;
      }
    }             

    hhl[Nz] = 0;
    for(int k=Nz;k>0;k--) 
      hhl[k-1] = hhl[k]+40.;
  }

  f2c_init_tenstream(fcomm,Nx,Ny,Nz, dx,dy, hhl, phi0, theta0,albedo);
  f2c_set_optical_properties(Nx,Ny,Nz, kabs, ksca, g);
  f2c_solve_tenstream(1.);
  get_result(Nx,Ny,Nz, edir,edn,eup,abso);

  f2c_destroy_tenstream();
  PetscFinalize();

  if (myid==0) {
    for(int j=0;j<Ny;j++){
      for(int k=0;k<Nz+1;k++)
        printf(" %f ",edir[k][j][0]);
      printf("\n");
      for(int k=0;k<Nz+1;k++)
        printf(" %f ",edn[k][j][0]);
      printf("\n");
      for(int k=0;k<Nz+1;k++)
        printf(" %f ",eup[k][j][0]);
      printf("\n      ");
      for(int k=0;k<Nz;k++)
        printf(" %f ",abso[k][j][0]);
      printf("\n");
      printf("\n");
    }
  }

  printf(" finished execution success!\n");

  MPI_Finalize();
  return(0);

}

