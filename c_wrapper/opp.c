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
#include <f2c_pprts.h>
#include "f2c_solver_ids.h"

#define CHKERR(int_status)                                                            \
  if (int_status) {                                                                   \
    fprintf(stderr, "error %d in %s ( %s::%d )", int_status, __func__, __FILE__, __LINE__); \
    return int_status; \
  }

// static char help[] = "This is the C wrapper interface to the pprts optprop environment.\n\n";

int main(int argc, char *argv[]) {

  if(argc>1) {
    fprintf(stderr, "Found cmd line args but did not expect any\n");
    for(int i=1; i < argc; ++i) {
      fprintf(stderr, "arg %d :: %s\n", i, argv[i]);
    }
  }

  const int mpi_comm  = 0;
  const int solver_id = SOLVER_ID_PPRTS_3_10;

  int ierr;
  int Ndir;
  int Ndiff;

  void * opp = NULL; // this is pointer to the tenstream optical properties object
  pprts_f2c_opp_init(mpi_comm, solver_id, &Ndir, &Ndiff, &opp, &ierr);
  CHKERR(ierr);

  fprintf(stderr, "LUT has %d direct and %d diffuse streams\n", Ndir, Ndiff);

  const float tauz        = 1;  // vertically integrated optical thickness
  const float w0          = .5; // single scatter albedo
  const float g           = .5; // asymmetry parameter
  const float aspect_zx   = 1;  // aspect ratio of a cell, i.e. dz / dx
  const float phi         = 0;  // solar azimuth angle
  const float theta       = 0;  // solar zenith angle
  const int lswitch_east  = 0;  // LUT's are only computed on angles [0,90] for azimuth
  const int lswitch_north = 0;  // set these values to switch internal couplings

  //
  // to see src and dst definitions of streams have a look in the tenstream source code
  // e.g. for 3/10 see `src/boxmc_3_10.inc`
  //
  { // direct 2 direct coeffs
    const int imode = 1;
    float dir2dir[Ndir][Ndir];
    pprts_f2c_opp_get_coeff(
        opp,
        tauz,
        w0,
        g,
        aspect_zx,
        phi,
        theta,
        imode,
        lswitch_east,
        lswitch_north,
        Ndir*Ndir,
        (float*)dir2dir,
        &ierr);
    CHKERR(ierr);

    fprintf(stdout, " * Direct to direct transport\n");
    fprintf(stdout, "      :: sources\n");
    for(int idst=0; idst<Ndir; ++idst) {
      fprintf(stdout, "dst %d ::", idst);
      for(int isrc=0; isrc<Ndir; ++isrc) {
        fprintf(stdout, " %f ", dir2dir[idst][isrc]);
      }
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
  }

  { // direct 2 diffuse coeffs
    const int imode = 2;
    float dir2diff[Ndiff][Ndir];
    pprts_f2c_opp_get_coeff(
        opp,
        tauz,
        w0,
        g,
        aspect_zx,
        phi,
        theta,
        imode,
        lswitch_east,
        lswitch_north,
        Ndiff*Ndir,
        (float*)dir2diff,
        &ierr);
    CHKERR(ierr);

    fprintf(stdout, " * Direct to diffuse transport\n");
    fprintf(stdout, "      :: sources\n");
    for(int idst=0; idst<Ndiff; ++idst) {
      fprintf(stdout, "dst %d ::", idst);
      for(int isrc=0; isrc<Ndir; ++isrc) {
        fprintf(stdout, " %f ", dir2diff[idst][isrc]);
      }
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
  }

  { // direct 2 diffuse coeffs
    const int imode = 3;
    float diff2diff[Ndiff][Ndiff];
    pprts_f2c_opp_get_coeff(
        opp,
        tauz,
        w0,
        g,
        aspect_zx,
        phi,
        theta,
        imode,
        lswitch_east,
        lswitch_north,
        Ndiff*Ndiff,
        (float*)diff2diff,
        &ierr);
    CHKERR(ierr);

    fprintf(stdout, " * Diffuse to diffuse transport\n");
    fprintf(stdout, "      :: sources\n");
    for(int idst=0; idst<Ndiff; ++idst) {
      fprintf(stdout, "dst %d ::", idst);
      for(int isrc=0; isrc<Ndiff; ++isrc) {
        fprintf(stdout, " %f ", diff2diff[idst][isrc]);
      }
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
  }
  pprts_f2c_opp_destroy(opp, &ierr);
  CHKERR(ierr);
}
