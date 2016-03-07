//-----Program CNU2D--------------------------------------------------
// (C) Copyright: Chungnam National University
//                Propeller & Cavitation Laboratory
//
//     Discrete Vortex/Source Method(DVM) Solution to
//       a Two Dimensional Steady Linear Airfoil Problem
//       without Using Explicit Kutta Condition
//
//     Force Calculation on 2-D Airfoil
//       with Discrete Vortex/Source Sept 76
//
//     CPP-Version 2.0: November 3, 1997  Gun-Do Kim
//                 2.1: Cleaned-up for Lecture Note, 2010, csl
//--------------------------------------------------------------------
#include <stdio.h>
#include <fstream.h>
#include <math.h>
#include "proplib.cpp"

const int kam2=17;  // kam2 = kam + 2; // num. of camber offsets.

//--- Function Prototypes ---
void indata();
void geom();
void setmat();
void force();

//--- Global Variables ---
int i, nt, ispace, *ipvt;
double alpha, pi, fc, tc, alphad, cosalf, sinalf, dtp;
double *xv, *zv, *xc, *zc, *xb, *sxb;
double *psq, *cctk, *cccm, *cam, *enx, *enz, *sig, *gam;
double **a, **alud;

// thickness and camber offset at 17-stations
static double per[] = {0.0 ,0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4,
                       0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 1.0};
static double thk[] = {0.0, 0.187, 0.2932, 0.4132, 0.5814, 0.8,
                       0.9274, 0.9904, 0.9924, 0.9306, 0.807, 0.622,
                       0.3754, 0.2286, 0.1496, 0.1018, 0.0};

//--- Define File Pointers ---
FILE  *fp_geo, *fp_gam, *fp_pre, *f_geo, *fp_mat;

void main()
{
    fp_geo = fopen("geom.tec", "w");
    fp_gam = fopen("gam.tec", "w");
    fp_pre = fopen("press.tec", "w");
    f_geo  = fopen("g.tec", "w");
    fp_mat = fopen("mat.out", "w");

    indata();
    geom();
    setmat();
    force();

    return;
}

//-------------------------------------------------------------------
//   indata(): read in the input
//-------------------------------------------------------------------
void indata()
{
    ifstream fin("foil.dat");
                    fin.ignore(80, '\n'); // Total num. of vortices(nt)
    fin >> nt; fin.ignore(80, '\n');
                    fin.ignore(80, '\n'); // Camber ratio (f/c)
    fin >> fc;  fin.ignore(80, '\n');
                    fin.ignore(80, '\n'); // Thick. ratio (t/c)
    fin >> tc;  fin.ignore(80, '\n'); 
                    fin.ignore(80, '\n'); // AoA in deg
    fin >> alphad;  fin.ignore(80, '\n');
                    fin.ignore(80, '\n'); // ispace=0: uniform sp
                    fin.ignore(80, '\n'); //       =1: cosine sp
                    fin.ignore(80, '\n'); //       =2: half-cosine
    fin >> ispace;  fin.ignore(80, '\n');

	printf("nt= %4d\n", nt);
    pi               = acos(-1.0);
    alpha            = alphad * pi / 180.0;
    cosalf           = cos(alpha);
    sinalf           = sin(alpha);

	//..allocate memory..
	xv  = new double[nt];   zv = new double[nt];
	xc  = new double[nt];   zc  = new double[nt];   xb = new double[nt+1];
	sxb = new double[nt+1]; psq = new double[kam2]; cam= new double[kam2];
	int mxspl = (kam2-1)*4;
	cctk= new double[mxspl];cccm= new double[mxspl];enx= new double[nt];
	enz = new double[nt];   sig=  new double[nt];   gam= new double[nt];

    return;
}

//-------------------------------------------------------------------
//   geom(): compute the location of vortex/source & control pts.
//-------------------------------------------------------------------
void geom()
{
    int i;
    double x_tilda;
    double *zup = new double[kam2];
    double *zlow= new double[kam2];

    for (i=0; i<kam2; i++) {
        // Parabolic meanline: y=4*x(1-x)
        cam[i]       = 4.0 * per[i] * (1.0 - per[i]);
        // printf("%8.6f\n", cam[i]);
        // Leading edge part is expanded by SQRT transformation
        // for interpolation of thickness 
        psq[i]       = sqrt(per[i]);
    }

    uglydk(kam2, 1, 1, psq, thk, 0.0, 0.0, cctk);
    uglydk(kam2, 1, 1, per, cam, 0.0, 0.0, cccm);

    if (ispace == 0) {  // uniform spacing
       for (i=0; i<nt; i++) {
           xv[i]        = (double(i) + 0.25) / double(nt);
           xb[i]        = double(i) / double(nt);
           sxb[i]       = sqrt(xb[i]);
           xc[i]        = xv[i] + 0.5 * ( 1.0 / double(nt) );
       }
           xb[nt]           = 1.0;
           sxb[nt]          = 1.0;
    }
    else if ( ispace == 1 ) {  // cosine spacing

       x_tilda  =  (pi / (double)nt);
           fprintf(f_geo, "variables = x, y\n");
           fprintf(f_geo, "zone i = %d\n", nt);
       for (i=0; i<nt; i++) {
           xb[i]      =  0.5 * ( 1.0 - cos(x_tilda*double(i)) );
           fprintf(f_geo,"%f  %f\n", xb[i], 1.0);
       }
           xb[nt]     =  1.0;
       for (i=0; i<nt; i++) {
           xv[i]      =  xb[i] + (( xb[i+1] - xb[i] ) * 0.25);
           xc[i]      =  xb[i] + (( xb[i+1] - xb[i] ) * 0.75);
       }
    }
    else { // half-cosine spacing

       x_tilda  =  (0.5 * pi / double(nt));
           fprintf(f_geo, "variables = x, y\n");
           fprintf(f_geo, "zone i = %d\n", nt);
       for (i=0; i<nt; i++) {
           xb[i]      = 1.0 -  cos(x_tilda*double(i));
       }
           xb[nt]     =  1.0;
       for (i=0; i<nt; i++) {
           xv[i]      =  xb[i] + (( xb[i+1] - xb[i] ) * 0.25);
           xc[i]      =  xb[i] + (( xb[i+1] - xb[i] ) * 0.75);
       }
       for (i=0; i<=nt; i++) {
            sxb[i]   = sqrt(xb[i]);
       }
    }

    //--- Interpolate vortex points on camber surface ---
    evaldk(kam2, nt, per, xv, zv, cccm);

    for (i=0; i<nt; i++ ) {
        zv[i]       *= fc;
    }

    evaldk(kam2, nt, per, xc, zc, cccm);

    for (i=0; i<nt; i++) {
        zc[i]       *= fc;
    }

    for (i=0; i<kam2; i++) {
        zup[i]       = cam[i] * fc  +  0.5 * thk[i] * tc;
        zlow[i]      = cam[i] * fc  -  0.5 * thk[i] * tc;
    }

    fprintf(fp_geo, "variables = x, z\n");
    fprintf(fp_geo, "zone i = %d\n", kam2);
    for (i=0; i<kam2; i++) {
        fprintf(fp_geo,"%f %f\n", per[i], zup[i]);
    }
    fprintf(fp_geo, "zone i = %d\n", kam2);
    for (i=0; i<kam2; i++) {
        fprintf(fp_geo,"%f %f\n", per[i], zlow[i]);
    }
    fprintf(fp_geo, "zone i = %d\n", kam2);
    for (i=0; i<kam2; i++) {
        fprintf(fp_geo,"%f %f\n", per[i], cam[i]*fc);
    }

    return;
}

//-------------------------------------------------------------------
//   setmat(): set coefficients of linear system
//-------------------------------------------------------------------
void setmat()
{
    int i, j;

    double slope, theta, dtpr, dx, dz;
    double units, wnits, unitv, wnitv, ludsign;
    double *tcn = new double[nt+1];
    double *dc  = new double[nt];
    double *d2c = new double[nt];
    double *zvi = new double[nt];
    double *ut  = new double[nt];
    double *wt  = new double[nt];

    //-- alloc memory for coefficient matrices ---
    ipvt= new int[nt];     
    a    = new double*[nt];
    for (i=0; i<nt; i++) {
        a[i] = new double[nt];
    }
    alud = new double*[nt];
    for (i=0; i<nt; i++) {
        alud[i] = new double[nt];
    }

    if ( ispace == 0 ) {
        evaldk(kam2, nt+1, psq, sxb, tcn, cctk);
    }
    else {
        evaldk(kam2, nt+1, per, xb, tcn, cctk);
    }

    for (i=0; i<nt; i++) {
        sig[i]       = tc * (tcn[i+1] - tcn[i]);
    }

    drivdk(kam2, nt, per, xc, dc, d2c, cccm);

    for (i=0; i<nt; i++) {
        slope        = fc * dc[i];
        theta        = atan(slope);
        enx[i]       = -sin(theta);
        enz[i]       = cos(theta);
    }

    dtp              = 0.5 / pi;
    for (i=0; i<nt; i++) {
        ut[i]        = 0.0;
        wt[i]        = 0.0;
        for (j=0; j<nt; j++) {
            dx       = xc[i] - xv[j];
            dz       = zc[i] - zv[j];
            dtpr    = dtp / (dx*dx + dz*dz);
            unitv    = dtpr * dz;
            wnitv    = -dtpr * dx;
            units    = dtpr * dx;
            wnits    = dtpr * dz;
            a[i][j]  = enx[i] * unitv + enz[i] * wnitv;
            ut[i]   += sig[j] * units;
            wt[i]   += sig[j] * wnits;
            fprintf(fp_mat, "%f   ", a[i][j]);
        }
        fprintf(fp_mat, " \n");
    }
    for (i=0; i<nt; i++) {
        gam[i] = -enx[i] * (cosalf+ut[i]) - enz[i] * (sinalf+wt[i]);
        fprintf(fp_mat, "%f   \n", gam[i]);
    }

    //--- LU decomposition of coeff. matrix ---
    ludcmp(a, alud, nt, ipvt, &ludsign);
    lubksb(alud, nt, ipvt, gam);

    fprintf(fp_gam, "variables = x, y\n");
    fprintf(fp_gam, "zone i = %d\n", nt);
    for (i=0; i<nt; i++) {
        fprintf(fp_gam, "%f %f\n", xc[i], gam[i]);
    }

    return;
}

//--------------------------------------------------------------------
//   force(): compute Lagally and Kutta-Joukowski forces
//--------------------------------------------------------------------
void force()
{
    int i, j;
    double tgam, dx, dz, bug, fxtt, fztt, fxlt, fzlt, force, cl,
          velup, velow, delnt;
    double *voncom=new double[nt];   double *vind = new double[nt];
    double *cpup  =new double[nt];   double *cplow= new double[nt];
    double *cpnega=new double[nt];   double *dgam = new double[nt];
    double *uts  = new double[nt];   double *wts  = new double[nt];
    double *uls  = new double[nt];   double *wls  = new double[nt];
    double *fxt  = new double[nt];   double *fzt  = new double[nt];
    double *fxl  = new double[nt];   double *fzl  = new double[nt];

    tgam             = 0.0;
    for (i=0; i<nt; i++) {
        tgam        += gam[i];
        uts[i]       = 0.0;
        wts[i]       = 0.0;
        uls[i]       = 0.0;
        wls[i]       = 0.0;
        for (j=0; j<nt; j++) {
            if (j == i) break;
            dx       = xv[i] - xv[j];
            dz       = zv[i] - zv[j];
            bug      = dtp / (dx*dx + dz*dz);
            uts[i]  += sig[j] * bug * dx;
            wts[i]  += sig[j] * bug * dz;
            uls[i]  += gam[j] * bug * dz;
            wls[i]  -= gam[j] * bug * dx;
        }
    }

    fxtt             = 0.0;
    fztt             = 0.0;
    fxlt             = 0.0;
    fzlt             = 0.0;
    for (i=0; i<nt; i++) {
        fxt[i]       = -sig[i] * (cosalf + uts[i] + uls[i]);
        fzt[i]       = -sig[i] * (sinalf + wts[i] + wls[i]);
        fxl[i]       = -gam[i] * (sinalf + wts[i] + wls[i]);
        fzl[i]       = gam[i] * (cosalf + uts[i] + uls[i]);
        fxtt        += fxt[i];
        fztt        += fzt[i];
        fxlt        += fxl[i];
        fzlt        += fzl[i];
    }
    force            = sqrt(pow((fxtt+fxlt),2) + pow((fztt+fzlt),2));
    cl               = 2.0 * force;

    printf("Lift Ceofficient (2-D VLM)   =  %f\n", cl);
	printf("Lift Ceofficient(2*pi*sin(alpha)) =  %f\n", 2*pi*sin(alpha));
	printf("Lift Ceofficient(4*pi*f/c) =  %f\n", 4*pi*fc);


    delnt            = 1.0 / double(nt);
    for (i=0; i<nt; i++) {
        voncom[i]    = cosalf * enz[i] - sinalf * enx[i];
        vind[i]      = gam[i] / (2.0 * double(delnt));
        velup        = voncom[i] + vind[i];
        velow        = voncom[i] - vind[i];
        cpup[i]      = 1.0 -  pow(velup,2);
        cplow[i]     = 1.0 -  pow(velow,2);
        cpnega[i]    = cpup[i] - cplow[i];
        dgam[i]      = gam[i] / double(delnt);
    }
    
	printf("Lift Ceofficient (2*Gam/Uc)   =  %f\n", 2*tgam);
    fprintf(fp_pre, "variables = x, z\n");
    fprintf(fp_pre, "zone i = %d\n", nt);
    for (i=0; i<nt; i++) {
        fprintf(fp_pre,"%f %f\n", xc[i], -cpup[i]);
    }
    fprintf(fp_pre, "zone i = %d\n", nt);
    for (i=0; i<nt; i++) {
        fprintf(fp_pre,"%f %f\n", xc[i], -cplow[i]);
    }
    return;
}
