#include <rayli_c_wrapper.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void VecSet(size_t N, float *arr, float val) {
    for(size_t i=0; i<N; i++) arr[i] = val;
}

int main() {

    int ierr;
    size_t Nphotons = 100000;
    size_t Nwedges = 1; // single wedge
    size_t Nfaces = 5;
    size_t Nverts = 6;
    int    cyclic = 1;
    float  *kabs, *ksca, *g, *albedo;
    float *flx_through_faces_edir, *flx_through_faces_ediff, *abso_in_cells;
    float sundir[] = {0,1/sqrt(2.),1/sqrt(2.)}; // vec towards the sun
    float diffuse_point_origin[] = {0.,0.,-1e30};

    kabs = malloc(Nwedges*sizeof(float)); VecSet(Nwedges, kabs, 0e-0);
    ksca = malloc(Nwedges*sizeof(float)); VecSet(Nwedges, ksca, 0e-0);
    g    = malloc(Nwedges*sizeof(float)); VecSet(Nwedges, g   , 1e-1);
    flx_through_faces_edir  = malloc(Nfaces*sizeof(float));
    flx_through_faces_ediff = malloc(2*Nfaces*sizeof(float));
    abso_in_cells           = malloc(Nwedges*sizeof(float));
    albedo                  = malloc(Nfaces*sizeof(float)); VecSet(Nfaces, albedo, -1);

    double vert_coords[]={ // 3 * Nverts
       0, 0, 1,       0, 0, 0,
       1, 0, 1,       1, 0, 0,
      .5, 0.866, 1,  .5, 0.866, 0
    };

    size_t verts_of_face[] = { // 4 * Nfaces
        0,4,2,-1,
        0,2,3,1,
        0,1,5,5,
        2,4,5,6,
        1,3,5,-1
    };

    size_t faces_of_wedges[] = { // 5 * Nwedges
      0,1,2,3,4
    };

    albedo[4] = .1;

    fprintf(stderr, "Huhu %s \n", rayli_version());

    ierr = rfft_wedge(Nphotons, Nwedges, Nfaces, Nverts, cyclic,
            verts_of_face, faces_of_wedges, vert_coords,
            kabs, ksca, g, albedo,
            sundir, diffuse_point_origin,
            flx_through_faces_edir,
            flx_through_faces_ediff,
            abso_in_cells);

    for(size_t f=0; f<Nfaces; f++) {
        fprintf(stderr, "on face %zu :: Edir %10g Edn %10g Eup %10g \n",
            f, flx_through_faces_edir[f], flx_through_faces_ediff[2*f], flx_through_faces_ediff[2*f+1]);
    }

    free(kabs);
    free(ksca);
    free(g);
    free(albedo);
    return ierr;
}


