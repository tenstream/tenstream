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

/*
 * Minimal toy reproducer for the master/slave optical-property distribution
 * path of the pprts C wrapper. It mirrors how libRadtran drives the solver:
 *
 *   - rank 0 (master) owns the *global* domain: it allocates and fills the
 *     full kabs/ksca/g/planck arrays and hands them to
 *     pprts_f2c_set_global_optical_properties(), which scatters them to the
 *     local MPI subdomains.
 *   - all other ranks (slaves) own *nothing*: they pass NULL for every field
 *     and only participate in the collective scatter/broadcast.
 *
 * This is exactly the call pattern that segfaults at np>=2 in a TenStream
 * built WITH_PETSC=OFF: the non-PETSc local_optprop() in
 * src/pprts.F90:set_global_optical_properties broadcasts the global arrays
 * via imp_bcast, but on the slave ranks those allocatable optionals are not
 * present -> dereferencing an absent argument -> crash. Run with:
 *
 *     mpirun -np 1 ./cwrapper_set_optprop      # works
 *     mpirun -np 2 ./cwrapper_set_optprop      # reproduces the slave crash
 */

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_PETSC
#include <petscsys.h>
#endif
#include <mpi.h>
#include <f2c_pprts.h>
#include "f2c_solver_ids.h"

#ifdef HAVE_PETSC
static char help[] = "Toy reproducer for pprts_f2c_set_global_optical_properties.\n\n";
#endif

static int collapseindex = 1;

/* rank 0: owns and distributes the global optical properties */
static int master(int fcomm) {
  int Nz = 16, Nx = 3, Ny = 4;
  int solver_id = SOLVER_ID_PPRTS_3_10;
  double dx = 100, dy = 100, dz = 50;
  float phi0 = 0, theta0 = 60;
  float albedo = 0.2;

  float *hhl    = (float *)malloc((Nz + 1) * sizeof(float));
  float *kabs   = (float *)malloc(Nz * Nx * Ny * sizeof(float));
  float *ksca   = (float *)malloc(Nz * Nx * Ny * sizeof(float));
  float *g      = (float *)malloc(Nz * Nx * Ny * sizeof(float));
  float *planck = (float *)malloc((Nz + 1) * Nx * Ny * sizeof(float));

  float *edir = (float *)malloc((Nz + 1) * Nx * Ny * sizeof(float));
  float *edn  = (float *)malloc((Nz + 1) * Nx * Ny * sizeof(float));
  float *eup  = (float *)malloc((Nz + 1) * Nx * Ny * sizeof(float));
  float *abso = (float *)malloc(Nz * Nx * Ny * sizeof(float));

  /* fill the global volume fields, layout is [Ny][Nx][Nz] (Nz fastest) */
  for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
      for (int k = 0; k < Nz; k++) {
        int ind = k + Nz * i + Nz * Nx * j;
        kabs[ind] = 1e-3;
        ksca[ind] = 1e-3;
        g[ind]    = 0.0;
      }

  /* planck lives on the Nz+1 level grid; keep it zero -> pure solar run */
  for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
      for (int k = 0; k < Nz + 1; k++)
        planck[k + (Nz + 1) * i + (Nz + 1) * Nx * j] = 0.0;

  /* monotone decreasing height of levels [m] */
  hhl[Nz] = 0;
  for (int k = Nz; k > 0; k--)
    hhl[k - 1] = hhl[k] + dz;

  pprts_f2c_init(fcomm, &solver_id, &Nz, &Nx, &Ny, &dx, &dy, hhl, &phi0, &theta0, &collapseindex);
  /* master passes the *global* arrays; the wrapper scatters them */
  pprts_f2c_set_global_optical_properties(Nz, Nx, Ny, &albedo, kabs, ksca, g, planck);
  pprts_f2c_solve(fcomm, 1.0);
  pprts_f2c_get_result(Nz, Nx, Ny, edn, eup, abso, edir);
  pprts_f2c_destroy(0);

  printf("edir[surface, i=0, j=0] = %f\n", edir[0]);
  printf("edn [surface, i=0, j=0] = %f\n", edn[0]);
  printf("eup [surface, i=0, j=0] = %f\n", eup[0]);

  free(hhl); free(kabs); free(ksca); free(g); free(planck);
  free(edir); free(edn); free(eup); free(abso);

  printf("master finished successfully\n");
  return 0;
}

/* all other ranks: own no global data, pass NULL everywhere */
static int slave(int fcomm) {
  int Nx, Ny, Nz;        /* filled in by pprts_f2c_init via broadcast */
  int solver_id = 0;
  double dx, dy;
  float phi0, theta0;

  float *albedo = NULL;
  float *kabs   = NULL;
  float *ksca   = NULL;
  float *g      = NULL;
  float *planck = NULL;
  float *hhl    = NULL;

  float *edir = NULL;
  float *edn  = NULL;
  float *eup  = NULL;
  float *abso = NULL;

  pprts_f2c_init(fcomm, &solver_id, &Nz, &Nx, &Ny, &dx, &dy, hhl, &phi0, &theta0, &collapseindex);
  /* slave passes NULL global arrays -> just joins the collective scatter */
  pprts_f2c_set_global_optical_properties(Nz, Nx, Ny, albedo, kabs, ksca, g, planck);
  pprts_f2c_solve(fcomm, 1.0);
  pprts_f2c_get_result(Nz, Nx, Ny, edn, eup, abso, edir);
  pprts_f2c_destroy(0);
  return 0;
}

int main(int argc, char *argv[]) {
  int numprocs, myid, fcomm;

  MPI_Init(&argc, &argv);
#ifdef HAVE_PETSC
  PetscInitialize(&argc, &argv, (char *)0, help);
  PetscInitializeFortran();
#endif

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
  fprintf(stderr, "rank %d of %d\n", myid, numprocs);

  if (myid == 0)
    master(fcomm);
  else
    slave(fcomm);

#ifdef HAVE_PETSC
  PetscFinalize();
#endif
  MPI_Finalize();
  return 0;
}
