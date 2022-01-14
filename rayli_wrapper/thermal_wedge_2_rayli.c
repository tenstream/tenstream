#include <rayli_c_wrapper.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void VecSet(size_t N, float *arr, float val) {
    for(size_t i=0; i<N; i++) arr[i] = val;
}

int main(int argc, char *argv[]) {

    int ierr;
    size_t Nthreads = 0;
    size_t Nphotons;
    size_t Nwedges = 2; // single wedge
    size_t Nfaces = 9;
    size_t Nverts = 8;
    int    cyclic = 1;
    float  *kabs, *ksca, *g, *albedo;
    double *flx_through_faces_ediff, *abso_in_cells;

    if(argc<=1) {
      Nphotons = 10000;
    } else {
      Nphotons = atoi(argv[1]);
    }

    kabs = malloc(Nwedges*sizeof(float)); VecSet(Nwedges, kabs, 1e+3);
    ksca = malloc(Nwedges*sizeof(float)); VecSet(Nwedges, ksca, 0e-0);
    g    = malloc(Nwedges*sizeof(float)); VecSet(Nwedges, g   , 1e-1);
    flx_through_faces_ediff = malloc(2*Nfaces*sizeof(double));
    abso_in_cells           = malloc(Nwedges*sizeof(double));
    albedo                  = malloc(Nfaces*sizeof(float)); VecSet(Nfaces, albedo, -1);

    double vert_coords[]={ // 3 * Nverts
       0, 0, 1,       0, 0, 0,
       1, 0, 1,       1, 0, 0,
       0, 1, 1,       0, 1, 0,
       1, 1, 1,       1, 1, 0
    };

    for (size_t i=0; i<3*Nverts; ++i) { vert_coords[i] *= 100; }

    //const float Bverts[]={ // Nverts
    //  100, 0,
    //  100, 0,
    //  100, 0,
    //  100, 0
    //};
    const float Bfaces[]={ // Nfaces
      100, 0, 0, 0, 10,
      100, 0, 0, 10
    };

    const size_t verts_of_face[] = { // 4 * Nfaces
        0,2,4,-1,
        0,2,3,1,
        0,1,5,4,
        2,4,5,3,
        1,3,5,-1,

        2,6,4,-1,
        4,5,7,6,
        2,6,7,3,
        3,7,5,-1
    };

    const size_t faces_of_wedges[] = { // 5 * Nwedges
      0,1,2,3,4,
      5,3,6,7,8
    };

    albedo[4] = .10;
    albedo[8] = .10;

    fprintf(stderr, "Huhu %s \n", rayli_version());

    ierr = rfft_wedge_thermal(Nthreads, Nphotons, Nwedges, Nfaces, Nverts, cyclic,
            verts_of_face, faces_of_wedges, vert_coords,
            kabs, ksca, g, albedo, Bfaces,
            flx_through_faces_ediff,
            abso_in_cells);
    if(ierr) {
      fprintf(stderr, "Rayli returned error %i\n", ierr);
      return ierr;
    }

    //for(size_t f=0; f<Nfaces; f++) {
    //    fprintf(stderr, "on face %zu :: Edn %10g Eup %10g \n",
    //        f, flx_through_faces_ediff[2*f], flx_through_faces_ediff[2*f+1]);
    //}
    size_t f;
    f = 0;
    fprintf(stderr, "wedge1 top face %zu :: Edn %10g Eup %10g \n", f, flx_through_faces_ediff[2*f], flx_through_faces_ediff[2*f+1]);
    f = 5;
    fprintf(stderr, "wedge2 top face %zu :: Edn %10g Eup %10g \n", f, flx_through_faces_ediff[2*f], flx_through_faces_ediff[2*f+1]);
    f = 4;
    fprintf(stderr, "wedge1 bot face %zu :: Edn %10g Eup %10g \n", f, flx_through_faces_ediff[2*f], flx_through_faces_ediff[2*f+1]);
    f = 8;
    fprintf(stderr, "wedge2 bot face %zu :: Edn %10g Eup %10g \n", f, flx_through_faces_ediff[2*f], flx_through_faces_ediff[2*f+1]);

    free(kabs);
    free(ksca);
    free(g);
    free(albedo);
    return ierr;
}


