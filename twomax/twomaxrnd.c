
/*----------------------------------------------------------------------
Full solution of the multi-layer two-stream equation (solar, thermal)
with maximum-random overlap assumption for partial cloudiness.
Author: Nina Črnivec, nina.crnivec@physik.uni-muenchen.de
----------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "twomaxrnd.h"

#define NFLUX 4   // number of diffuse fluxes defined at each level (Edn_c, Edn_f, Eup_c, Eup_f)
#define GAUSS_SINGULAR  -20

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef MAX
#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#endif
#ifndef MIN
#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

/*
OUTPUT:
Edir = direct irradiance;
Edn  = diffuse downward irradiance;
Eup  = diffuse upward irradiance;
*/

// Effective thermal emission of a layer with planck emission values defined at the boundaries
double B_eff_mu(double B_far, double B_near, double tau, double T, double mu) {
    if (tau/mu<1.e-8) {
        return (B_far+B_near)*.5;
    } else {
        return (-B_near + B_far * T)/(-1. + T) + ((B_far - B_near) *mu)/tau;
    }
    return -1;
}
double B_eff(double B_far, double B_near, double tau) {
    const int Nmu=20;
    const double dmu = 1./Nmu;
    double Beff=0;
    for (int imu = 0; imu<Nmu; ++imu) {
        double mu = .5*dmu + imu*dmu;
        Beff += B_eff_mu(B_far, B_near, tau, exp(-tau/mu), mu) * mu * dmu;
    }
    return Beff*2;
}

/***********************************************************************************/
/* Function: solve_gauss                                                  @31_30i@ */
/* Description:                                                                    */
/*  Solve a system of n linear equations, A*x = b, using the Gauss algorithm       */
/*  The pivot element is determined using 'relative column maximum strategy'.      */
/*  For a description of the algorithm see H.R.Schwarz: "Numerische Mathematik",   */
/*  pg. 21. Memory for the result vector res is allocated automatically.           */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double **A:            Matrix[n x n]  (see above).                             */
/*  double *b:             Vector[n] (see above).                                  */
/*  double n:              Number of equations.                                    */
/*  double **res:          Pointer to the result vector[n]; if no unique solution  */
/*                         exists, *res will be set to NULL.                       */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if no unique solution.                                          */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

static int solve_gauss (double **A, double *b, int n, double *res)
{
    int i=0, j=0, k=0, p=0;
    double div=0.0;
    double *tmp1=NULL;
    double tmp2=0.0;
    double sum=0.0;
    double max = (0.0 - DBL_MAX);

    for (k=0; k<n; k++)  {
        /* get pivot element (line p) */
        max = (0.0 - DBL_MAX);
        p = 0;

        for (i=k; i<n; i++)  {

            sum = 0.0;
            for (j=k; j<n; j++)
                sum += fabs(A[i][j]);
            if (sum == 0.0)
                return GAUSS_SINGULAR;

            if ( (tmp2 = fabs (A[i][k]) / sum) > max )  {
                max = tmp2;
                p = i;
            }
        }

        /* exchange lines k and p */
        tmp1 = A[k];
        A[k] = A[p];
        A[p] = tmp1;

        tmp2 = b[k];
        b[k] = b[p];
        b[p] = tmp2;


        if ( (div = A[k][k]) == 0)   /* no definite solution */
            return GAUSS_SINGULAR;

        for (i=k; i<n; i++)
            A[k][i] /= div;

        b[k] /= div;

        for (i=k+1; i<n; i++)  {
            div = A[i][k];
            for (j=k; j<n; j++)
                A[i][j] = A[i][j] - A[k][j] * div;
            b[i] = b[i] - b[k] * div;
        }
    }

    for (i=0; i<n; ++i) res[i] = 0;

    for (i=n-1; i>=0; i--)  {
        res[i] += b[i];
        for (k=i+1; k<n; k++)
            res[i] -= A[i][k] * res[k];
    }

    return 0;
}


//======================
// FUNCTION calcp1p2p3p4
//======================
// INPUT:  nlev, cf;
// OUTPUT: ar_p1, ar_p2, ar_p3, ar_p4;

/*
 * From Zdunkowski (pages 180-183):
 * p1 = (1 - max(cf[i], cf[i-1])) / (1 - cf[i-1])
 * p2 = (1 - max(cf[i], cf[i+1]) / (1 - cf[i+1])
 * p3 = min(cf[i], cf[i-1]) / cf[i-1]
 * p4 = min(cf[i], cf[i+1]) / cf[i+1]

 * Remark from Zdunkowski:
 * It is noteworthy that a particular coefficient pj (j=1,2,3,4) is set equal to 1
 * if an undetermined expression 0/0 occurs.
 * This follows from physical reasoning or from applying l'Hopital's rule.
 */
static int calcp1p2p3p4 (int nlev, double *cf, double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4)
{
    double max_p1; // = max(cf[ilyr],cf[ilyr-1]); appears in the expression for p1;
    double min_p3; // = min(cf[ilyr],cf[ilyr-1]); appears in the expression for p3;
    double max_p2; // = max(cf[ilyr],cf[ilyr+1]); appears in the expression for p2;
    double min_p4; // = min(cf[ilyr],cf[ilyr+1]); appears in the expression for p4;

    const int nlyr=nlev-1;

    /*
     * FORMULAS:
     * ar_p1[ilyr] = (1.0 - max(cf[ilyr], cf[ilyr-1])) / (1.0 - cf[ilyr-1]);
     * ar_p2[ilyr] = (1.0 - max(cf[ilyr], cf[ilyr+1]) / (1.0 - cf[ilyr+1]);
     * ar_p3[ilyr] = min(cf[ilyr], cf[ilyr-1]) / cf[ilyr-1];
     * ar_p4[ilyr] = min(cf[ilyr], cf[ilyr+1]) / cf[ilyr+1];
     */

    // Calculate vertical profiles of p1, p3:
    // Special case: ilyr=0;
    ar_p1[0] = 1.0;
    ar_p3[0] = 1.0;
    for(int ilyr=0;ilyr<nlyr;ilyr++){
            // Find max_p1:
            if(cf[ilyr] > cf[ilyr-1]){
                max_p1 = cf[ilyr];
            }else{
                max_p1 = cf[ilyr-1];
            }//e-if

            // Find_min_p3:
            if(cf[ilyr] < cf[ilyr-1]){
                min_p3 = cf[ilyr];
            }else{
                min_p3 = cf[ilyr-1];
            }//e-if

            if(cf[ilyr-1] == 1.0) ar_p1[ilyr] = 1.0;
            else ar_p1[ilyr] = (1.0 - max_p1) / (1.0 - cf[ilyr-1]);

            if(cf[ilyr-1] == 0.0) ar_p3[ilyr] = 1.0;
            else ar_p3[ilyr] = min_p3 / cf[ilyr-1];
    }//e-for over ilyr;

    // Calculate vertical profiles of p2 and p4:
    // Special case: ilyr = nlyr-1 = nlev-2;
    for(int ilyr=0;ilyr<nlyr-1;ilyr++){

            // Find max_p2:
            if(cf[ilyr] > cf[ilyr+1]){
                max_p2 = cf[ilyr];
            }else{
                max_p2 = cf[ilyr+1];
            }//e-if

            // Find_min_p4:
            if(cf[ilyr] < cf[ilyr+1]){
                min_p4 = cf[ilyr];
            }else{
                min_p4 = cf[ilyr+1];
            }//e-if

            if(cf[ilyr+1] == 1.0) ar_p2[ilyr] = 1.0;
            else ar_p2[ilyr] = (1.0 - max_p2) / (1.0 - cf[ilyr+1]);

            if(cf[ilyr+1] == 0.0) ar_p4[ilyr] = 1.0;
            else ar_p4[ilyr] = min_p4 / cf[ilyr+1];
    }//e-for
    ar_p2[nlyr-1] = 1.0;
    ar_p4[nlyr-1] = 1.0;

    return 0;
}//e-calcp1p2p3p4



//===========================
// FUNCTION: eddington_coeffc
//===========================

/*
   Calculate Eddington coefficients a11, a12, a13, a23, and a33
   from layer optical thickness dtau, asymmetry parameter g,
   single scattering albedo omega0, and cosine of solar zenith angle mu0.
   */

static void eddington_coeffc_zdun (double dtau, double g, double omega0, double mu0,
        double *a11, double *a12, double *a13, double *a23, double *a33)
{
    const double eps_resonance = 1e-8;
    dtau = MIN(200, MAX(1e-16, dtau));
    g    = MAX(1e-6, g);
    omega0 = MIN(1.-eps_resonance, MAX(1e-16,omega0));

    double mubar = .5;
    double bbar  = 3./8.*(1.-g);
    double b_minus_mu0 = .5 - .75 * g *mu0;

    double alpha_1 = ( 1.- omega0*(1.-bbar) ) / mubar;
    double alpha_2 = omega0*bbar/mubar;

    double bscr = 0.5 - 0.375 * g;
    alpha_1 = 2. * ( 1. - omega0 * ( 1. - bscr ) ) - 0.25;
    alpha_2 = 2. * omega0 * bscr - 0.25;

    double lambda = sqrt(alpha_1*alpha_1 - alpha_2*alpha_2);

    double e1 = exp( MIN(DBL_MAX_EXP, lambda*dtau));
    double e2 = exp(-MIN(DBL_MAX_EXP, lambda*dtau));

    double alpha1_m_lambda = alpha_1-lambda ;
    double alpha1_p_lambda = alpha_1+lambda ;
    if((alpha1_m_lambda>-1e-16) && (alpha1_m_lambda<1e-16))
        alpha1_m_lambda = copysign(1e-16, alpha1_m_lambda);

    if((alpha1_p_lambda>-1e-16) && (alpha1_p_lambda<1e-16))
        alpha1_p_lambda = copysign(1e-16, alpha1_p_lambda);

    double A = 1. / ( alpha_2/alpha1_m_lambda*e1 - alpha_2/alpha1_p_lambda * e2 );

    double beta11  =  A * alpha_2/alpha1_m_lambda;
    double beta21  = -A * alpha_2/alpha1_p_lambda;
    double beta12  = -A * e2;
    double beta22  =  A * e1;

    double gamma12 = alpha_2/alpha1_p_lambda * e1;
    double gamma22 = alpha_2/alpha1_m_lambda * e2;

    *a11 = beta11 + beta21;
    *a12 = beta12 + beta22;

    *a11 = MAX(0.,  MIN(1., *a11) );
    *a12 = MAX(0.,  MIN(1., *a12) );

    if(mu0>1e-16) {
        *a33     = exp ( - MIN(DBL_MAX_EXP, dtau / mu0 ));

        double alpha_3 = -omega0 * b_minus_mu0;
        double alpha_4 =  omega0 * (1.-b_minus_mu0);
        double den = (1./mu0/mu0) - lambda*lambda ;
        if( fabs(den)<eps_resonance ) {
            if(mu0<.5) {
                den = 1./ (mu0*mu0 - eps_resonance)  - lambda*lambda;
            }else{
                den = 1./ (mu0*mu0 + eps_resonance)  - lambda*lambda;
            }
        }

        double alpha_5 = ( (alpha_1-1./mu0)*alpha_3 - alpha_2*alpha_4 ) / den;
        double alpha_6 = ( alpha_2*alpha_3 - (alpha_1+1./mu0)*alpha_4 ) / den;

        double beta13  = -beta11*alpha_5 * *a33 - beta12*alpha_6;
        double beta23  = -beta21*alpha_5 * *a33 - beta22*alpha_6;

        *a13 = beta13         + beta23         + alpha_5;
        *a23 = beta13*gamma12 + beta23*gamma22 + alpha_6 * *a33;

        //*a13 = *a13 / mu0 !Fabian: Roberts coefficients a13 expect S to be
        //*a23 = *a23 / mu0 !        irradiance on tilted plane... we use irradiance on z-plane

        *a13 = MAX(0., *a13);
        *a23 = MAX(0., *a23);
    }else{
        *a33=0;
        *a13=0;
        *a23=0;
    }
}//e-eddington_coeffc


//=======================
// FUNCTION: buildMatrixA
//=======================
// INPUT:  nlev, Ag, ar_a11_c, ar_a11_f, ar_a12_c, ar_a12_f, ar_p1, ar_p2, ar_p3, ar_p4;
// OUTPUT: matrixA;

static int buildMatrixA (int nlev, double Ag,
        double *ar_a11_c, double *ar_a11_f,
        double *ar_a12_c, double *ar_a12_f,
        double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4,
        double **matrixA)
{
    int i;     // position in levels (e.g. for nlev=21, i in range: 0 - 20)
    int iRow;  // row position in matrix A (e.g. for nlev=21, iRow in range: 0 - 83)

    int nextColEup; // index of column for next Eup data
    int nextColEdw; // index of column for next Edn data; nextColEdn = nextColEup - 4

    // Set values for initial four rows of matrix A:
    // First row:
    matrixA[0][2] = ar_a12_f[0]*ar_p1[0];
    matrixA[0][3] = ar_a12_f[0]*(1.0-ar_p3[0]);
    matrixA[0][4] = ar_a11_f[0]*ar_p2[0];
    matrixA[0][5] = ar_a11_f[0]*(1.0-ar_p4[0]);

    // Second row:
    matrixA[1][2] = ar_a12_c[0]*(1.0-ar_p1[0]);
    matrixA[1][3] = ar_a12_c[0]*ar_p3[0];
    matrixA[1][4] = ar_a11_c[0]*(1.0-ar_p2[0]);
    matrixA[1][5] = ar_a11_c[0]*ar_p4[0];
    // Third row is already zero; (needs to be zero due to upper boundary condition);
    // Forth row is already zero; (needs to be zero due to upper boundary condition);

    nextColEup = 6; // index of column for Eup(level1)
    nextColEdw = 2; // index of column for Edn(level1)

    for(i=1;i<nlev;i++){
        iRow = NFLUX*i;

        if(i == (nlev-1)){ // lower boundary condition for the forth and third row from bottom up of matrix A;
            matrixA[iRow][nextColEup] = Ag;
            matrixA[iRow+1][nextColEup+1] = Ag;
        }else{
            matrixA[iRow][nextColEup]   = ar_a12_f[i]*ar_p1[i];
            matrixA[iRow][nextColEup+1] = ar_a12_f[i]*(1.0-ar_p3[i]);
            matrixA[iRow][nextColEup+2] = ar_a11_f[i]*ar_p2[i];
            matrixA[iRow][nextColEup+3] = ar_a11_f[i]*(1.0-ar_p4[i]);

            matrixA[iRow+1][nextColEup]   = ar_a12_c[i]*(1.0-ar_p1[i]);
            matrixA[iRow+1][nextColEup+1] = ar_a12_c[i]*ar_p3[i];
            matrixA[iRow+1][nextColEup+2] = ar_a11_c[i]*(1.0-ar_p2[i]);
            matrixA[iRow+1][nextColEup+3] = ar_a11_c[i]*ar_p4[i];
        }//e-if
        nextColEup += NFLUX;

        matrixA[iRow+2][nextColEdw]   = ar_a11_f[i-1]*ar_p1[i-1];
        matrixA[iRow+2][nextColEdw+1] = ar_a11_f[i-1]*(1.0-ar_p3[i-1]);
        matrixA[iRow+2][nextColEdw+2] = ar_a12_f[i-1]*ar_p2[i-1];
        matrixA[iRow+2][nextColEdw+3] = ar_a12_f[i-1]*(1.0-ar_p4[i-1]);

        matrixA[iRow+3][nextColEdw]   = ar_a11_c[i-1]*(1.0-ar_p1[i-1]);
        matrixA[iRow+3][nextColEdw+1] = ar_a11_c[i-1]*ar_p3[i-1];
        matrixA[iRow+3][nextColEdw+2] = ar_a12_c[i-1]*(1.0-ar_p2[i-1]);
        matrixA[iRow+3][nextColEdw+3] = ar_a12_c[i-1]*ar_p4[i-1];

        nextColEdw += NFLUX;
    }//e-for

    return 0;
}//e-buildMatrixA


//=========================
// FUNCTION buildVectorBsol
//=========================
// INPUT:  nlev, Ag, mu0, ar_a13_c, ar_a13_f, ar_a23_c, ar_a23_f, ar_S_c, ar_S_f, ar_p1, ar_p3;
// OUTPUT: vectB;

static int buildVectorBsol (int nlev, double Ag, double mu0,
        double *ar_a13_c, double *ar_a13_f,
        double *ar_a23_c, double *ar_a23_f,
        double *ar_S_c, double *ar_S_f,
        double *ar_p1, double *ar_p3,
        double *vectB)
{
    int i;  // position in levels
    int j;  // position in vector B

    // Set initial four values:
    vectB[0] = ar_a13_f[0]*ar_S_f[0];
    vectB[1] = ar_a13_c[0]*ar_S_c[0];
    vectB[2] = 0.0; // upper boundary condition
    vectB[3] = 0.0; // upper boundary condition

    for(i=1; i<(nlev-1);i++){
        j = NFLUX*i;
        vectB[j]   = ar_a13_f[i]*(ar_p1[i]*ar_S_f[i] + (1.0-ar_p3[i])*ar_S_c[i]);
        vectB[j+1] = ar_a13_c[i]*((1.0-ar_p1[i])*ar_S_f[i] + ar_p3[i]*ar_S_c[i]);
        vectB[j+2] = ar_a23_f[i-1]*(ar_p1[i-1]*ar_S_f[i-1] + (1.0-ar_p3[i-1])*ar_S_c[i-1]);
        vectB[j+3] = ar_a23_c[i-1]*((1.0-ar_p1[i-1])*ar_S_f[i-1] + ar_p3[i-1]*ar_S_c[i-1]);
    }//e-for

    // Treat last four values seperately:
    j = NFLUX*(nlev-1);
    vectB[j]   = Ag*mu0*ar_S_f[nlev-1]; // lower boundary condition
    vectB[j+1] = Ag*mu0*ar_S_c[nlev-1]; // lower boundary condition
    vectB[j+2] = ar_a23_f[nlev-2]*(ar_p1[nlev-2]*ar_S_f[nlev-2] + (1.0-ar_p3[nlev-2])*ar_S_c[nlev-2]);
    vectB[j+3] = ar_a23_c[nlev-2]*((1.0-ar_p1[nlev-2])*ar_S_f[nlev-2] + ar_p3[nlev-2]*ar_S_c[nlev-2]);

    return 0;
}//e-buildVectorBsol



//=========================
// FUNCTION buildVectorBthe
//=========================
// INPUT:  nlev, Ag, Bg, ar_theComp1_c, ar_theComp1_f, ar_theComp2_c, ar_theComp2_f, cf;
// OUTPUT: bb_the;

// thermal component 1 = upward emission
// thermal component 2 = downward emission

static int buildVectorBthe (int nlev, double Ag, double Bg,
        double *ar_theComp1_c, double *ar_theComp1_f,
        double *ar_theComp2_c, double *ar_theComp2_f,
        double *cf,
        double *vectB)
{
    int i;  // position in levels
    int j;  // position in vector B

    // Set initial four values:
    vectB[0] = (1.0-cf[0])*ar_theComp1_f[0];
    vectB[1] = cf[0]*ar_theComp1_c[0];
    vectB[2] = 0.0; // upper boundary condition
    vectB[3] = 0.0; // upper boundary condition

    for(i=1; i<(nlev-1); i++){
        j = NFLUX*i;
        vectB[j]   = (1.0-cf[i])*ar_theComp1_f[i];
        vectB[j+1] = cf[i]*ar_theComp1_c[i];
        vectB[j+2] = (1.0-cf[i-1])*ar_theComp2_f[i-1];
        vectB[j+3] = cf[i-1]*ar_theComp2_c[i-1];
    }//e-for

    // Treat last four values separately:
    // i=nlev-1; // bottom (ground) level;
    j = NFLUX*(nlev-1);
    vectB[j]   = (1.0-cf[nlev-2])*(1.0-Ag)*M_PI*Bg; // lower boundary condition
    vectB[j+1] = cf[nlev-2]*(1.0-Ag)*M_PI*Bg;       // lower boundary condition
    vectB[j+2] = (1.0-cf[nlev-2])*ar_theComp2_f[nlev-2];
    vectB[j+3] = cf[nlev-2]*ar_theComp2_c[nlev-2];

    return 0;
}//e-buildVectorBthe


int twostream_maxrand (double *dtau_org_c, double *omega0_org_c, double *g_org_c, //"_c" for cloudy region
        double *dtau_org_f, double *omega0_org_f, double *g_org_f, //"_f" for cloud-free region
        double *cf, int nlev, double S0, double mu0, double Ag,
        double Bg, double *B, int delta, int flagSolar, int flagThermal,
        double *Edir, double *Edn, double *Eup)
{
    const int nlyr=nlev-1;
    int iStatus=0;

    /*
       Eddington coefficients:
       a11 = transmission coefficient for diffuse radiation;
       a12 = reflection coefficient for diffuse radiation;
       a13 = reflection coefficient for the primary scattered parallel solar radiation;
       a23 = transmission coefficient for the primary scattered parallel solar radiation;
       a33 = transmission coefficient for the direct parallel solar radiation;
       */
    double a11_c[nlyr];
    double a12_c[nlyr];
    double a13_c[nlyr];
    double a23_c[nlyr];
    double a33_c[nlyr];

    double a11_f[nlyr];
    double a12_f[nlyr];
    double a13_f[nlyr];
    double a23_f[nlyr];
    double a33_f[nlyr];

    // Components of vector B in the thermal spectral range:
    double theComp1_c[nlev];
    double theComp1_f[nlev];
    double theComp2_c[nlev];
    double theComp2_f[nlev];

    // Parameters related to cloud cover of two adjacent layers:
    // in Zdunkowski (pages: 180-183) denoted as: b1,b2,b3,b4;
    // here: p1, p2, p3, p4;
    double p1[nlyr];
    double p2[nlyr];
    double p3[nlyr];
    double p4[nlyr];

    double bb_sol[NFLUX*nlev];
    double bb_the[NFLUX*nlev];
    double bb    [NFLUX*nlev];
    double xx    [NFLUX*nlev]; // result vector

    double S_c[nlev];
    double S_f[nlev];

    double Edir_c[nlev];
    double Edir_f[nlev];

    double Eup_c[nlev];
    double Eup_f[nlev];

    double Edn_c[nlev];
    double Edn_f[nlev];

    double **AA=NULL;

    for(int k=0;k<nlev-1;++k){
        int ierr=0;

        if(g_org_f[k]>1 || g_org_f[k]<0)
        {printf("Found bad value in clear sky layer %d for assym param g: %f \n", k, g_org_f[k]); ierr++; }

        if(g_org_c[k]>1 || g_org_c[k]<0)
        {printf("Found bad value in cloudy    layer %d for assym param g: %f \n", k, g_org_c[k]); ierr++; }

        if(omega0_org_f[k]>1 || omega0_org_f[k]<0)
        {printf("Found bad value in clear sky layer %d for w0: %f \n", k, omega0_org_f[k]); ierr++; }

        if(omega0_org_c[k]>1 || omega0_org_c[k]<0)
        {printf("Found bad value in cloudy    layer %d for w0: %f \n", k, omega0_org_c[k]); ierr++; }

        if(ierr!=0) return ierr;
    }
    AA=malloc(NFLUX*nlev*sizeof(double*));

    for(int ilev=0;ilev<NFLUX*nlev;ilev++){
        if((AA[ilev]=calloc(NFLUX*nlev, sizeof(double)))==NULL){
            fprintf (stderr, "Error allocating memory for AA[%d]\n", ilev);
            return -1;
        }//e-if
    }//e-for

    // At the moment it is only possible to calculate solar OR thermal RT, but not both simultaneously;
    if ((flagSolar && flagThermal) || (!flagSolar && !flagThermal)){
        fprintf (stderr, "Error - invalid input parameters - use flagSolar and flagThermal alternatingly\n");
        fprintf (stderr, "flagSolar = %d , flagThermal = %d\n", flagSolar, flagThermal);
        return -2;
    }//e-if

    // Calculate vertical profiles of p1, p2, p3 and p4 from vertical profile of cloud fraction:
    iStatus = calcp1p2p3p4(nlev, cf, p1, p2, p3, p4);

    if(iStatus != 0){
        fprintf (stderr, "Error calculating vertical profiles of p1, p2, p3 and p4 from cloud cover; ERROR=%d \n", iStatus);
        return -3;
    }//e-if


    // Calculate vertical profiles of Eddington coefficients
    // and vertical profiles of thermal components:
    for(int k=0;k<nlyr;k++){
        eddington_coeffc_zdun (
                dtau_org_c[k],
                g_org_c[k],
                omega0_org_c[k],
                mu0,
                &a11_c[k],
                &a12_c[k],
                &a13_c[k],
                &a23_c[k],
                &a33_c[k]);

        eddington_coeffc_zdun (dtau_org_f[k],
                g_org_f[k],
                omega0_org_f[k],
                mu0,
                &a11_f[k],
                &a12_f[k],
                &a13_f[k],
                &a23_f[k],
                &a33_f[k]);

        if(flagThermal){
            theComp1_c[k] = B_eff(B[k+1], B[k], dtau_org_c[k])*M_PI*(1.-a11_c[k]);
            theComp2_c[k] = B_eff(B[k], B[k+1], dtau_org_c[k])*M_PI*(1.-a11_c[k]);
            theComp1_f[k] = B_eff(B[k+1], B[k], dtau_org_f[k])*M_PI*(1.-a11_f[k]);
            theComp2_f[k] = B_eff(B[k], B[k+1], dtau_org_f[k])*M_PI*(1.-a11_f[k]);
        }//e-if
    }//e-for

    // Initialize vectors S_c and S_f:
    // S[0]= S0 = S_c[0] + S_f[0] = cf[0]*S0 + (1.0-cf[0])*S0;
    if(flagSolar){
        S_c[0]=cf[0]*S0;
        S_f[0]=(1.0-cf[0])*S0;

        for(int ilev=1;ilev<nlev;ilev++){
            S_c[ilev]=a33_c[ilev-1]*((1.0-p1[ilev-1])*S_f[ilev-1] + p3[ilev-1]*S_c[ilev-1]);
            S_f[ilev]=a33_f[ilev-1]*(p1[ilev-1]*S_f[ilev-1] + (1.0-p3[ilev-1])*S_c[ilev-1]);
        }//e-for
    }//e-if

    // Equation system has the following form: xx = AA*xx + bb;
    // Build matrix AA:
    iStatus = buildMatrixA (nlev, Ag, a11_c, a11_f, a12_c, a12_f,
            p1, p2, p3, p4, AA);
    if(iStatus != 0){
        fprintf (stderr, "buildMatrixA ERROR=%d \n", iStatus);
        return -4;
    }//e-if

    // Make matrix A2 = AA - IdentityMatrix and save it to the same matrix AA,
    // since the equation system in the form Ax=b must be passed to solve_gauss;
    // xx = AA*xx + bb;
    // (AA-II)*xx = -bb;
    for(int iRow=0; iRow < NFLUX*nlev; iRow++) AA[iRow][iRow] = -1.0;

    if(iStatus != 0){
        fprintf (stderr, "makeMatrixA2 ERROR=%d \n", iStatus);
        return -5;
    }//e-if

    // Build vector bb_sol:
    if(flagSolar){
        iStatus = buildVectorBsol (nlev, Ag, mu0, a13_c, a13_f, a23_c, a23_f, S_c, S_f, p1, p3, bb_sol);

        if(iStatus != 0){
            fprintf (stderr, "buildVectorBsol ERROR=%d \n", iStatus);
            return -6;
        }//e-if
    }//e-if

    // Build vector bb_the:
    if(flagThermal){
        iStatus = buildVectorBthe (nlev, Ag, Bg, theComp1_c, theComp1_f, theComp2_c, theComp2_f, cf, bb_the);

        if(iStatus != 0){
            fprintf (stderr, "buildVectorBthe ERROR=%d \n", iStatus);
            return -7;
        }//e-if
    }//e-if

    // Assign bb_sol or bb_the to the final vector bb:
    for(int k=0; k<NFLUX*nlev; k++) bb[k] = 0;
    if(flagSolar)   { for(int k=0; k<NFLUX*nlev; k++) bb[k] += bb_sol[k]; }
    if(flagThermal) { for(int k=0; k<NFLUX*nlev; k++) bb[k] += bb_the[k]; }

    // Make vector -bb and save it to the same vector bb:
    for(int iRow=0; iRow < NFLUX*nlev; iRow++) bb[iRow] = -bb[iRow];

    if(iStatus != 0){
        fprintf (stderr, "makeVectorMinusB ERROR=%d \n", iStatus);
        return -8;
    }//e-if

    // Solve AA*xx=bb to obtain xx:
    // AA is a 11-diagonal matrix; calculation could be sped up with a 11
    iStatus = solve_gauss (AA, bb, NFLUX*nlev, xx);

    if(iStatus != 0){
        fprintf (stderr, "Error %d solving equation system\n", iStatus);
        return -9;
    }//e-if

    for(int ilev=0;ilev<nlev;ilev++){
        Eup_f[ilev] = xx[NFLUX*ilev];
        Eup_c[ilev] = xx[NFLUX*ilev+1];
        Edn_f[ilev] = xx[NFLUX*ilev+2];
        Edn_c[ilev] = xx[NFLUX*ilev+3];
    }//e-for

    for(int ilev=0;ilev<nlev;ilev++){
        Edir_c[ilev] = S_c[ilev]*mu0;
        Edir_f[ilev] = S_f[ilev]*mu0;
    }//e-for

    // Sum up the irradiances for cloudy and cloud-free regions to obtain the final result:
    for(int ilev=0;ilev<nlev;ilev++){
        Edir[ilev] = Edir_c[ilev] + Edir_f[ilev];
        Eup[ilev]  = Eup_c[ilev] + Eup_f[ilev];
        Edn[ilev]  = Edn_c[ilev] + Edn_f[ilev];
    }//e-for

    for(int k=0; k < NFLUX*nlev; k++) free(AA[k]);
    free(AA);

    return iStatus;
}//e-twostream_maxrand
