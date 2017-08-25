void f2c_tenstream_rrtmg(int fcomm,
    int *Nz, int *Nx,int *Ny,
    double *dx, double *dy,
    double *phi0, double *theta0,
    double *albedo_thermal, double* albedo_solar,
    char *atm_filename,
    int *lthermal, int *lsolar,
    int *Nz_merged, double **edir,
    double **edn, double **eup, double **abso,
    double *d_plev, double *d_tlev,
    double *d_lwc, double *d_reliq,
    double *d_iwc, double *reice,
    int *nprocx,
    int *nxproc,
    int *nprocy,
    int *nyproc
    );

void f2c_destroy_tenstream_rrtmg(int *lfinalizepetsc);
