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




///////////////   Program CNU2D Case 8 make by SeokHo ////////////////
//--------------------------------------------------------------------
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "proplib1.h"

#define PI 3.1415926535

double realsin( double );
double realcos( double );

const int kam2=17;  // kam2 = kam + 2; // num. of camber offsets.
using namespace std;
//--- Function Prototypes ---
void indata();
void geom();
void setmat();
void force();

//--- Global Variables ---
int i, nt, ispace, *ipvt, *subipvt;

double leftright, updown;

double alpha, subalpha, pi, fc, tc, alphad, subalphad, cosalf, sinalf, subcosalf, subsinalf, dtp, subdtp;
double *xv, *zv, *xc, *zc, *xb, *sxb;
double *psq, *cctk, *cccm, *cam, *enx, *enz, *sig, *gam;

double *subxv, *subzv, *subxc, *subzc, *subxb, *subsxb;
double *subpsq, *subcctk, *subcccm, *subcam, *subenx, *subenz, *subsig, *subgam;

double *pxv, *pzv, *pxc, *pzc, *pxb, *psxb;
double *pcam,*pgam;

double *psubxv, *psubzv, *psubxc, *psubzc, *psubxb, *psubsxb;
double *psubcam, *psubgam;

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
    fp_geo = fopen("geom.dat", "w");
    fp_gam = fopen("gam.dat", "w");
    fp_pre = fopen("press.dat", "w");
    f_geo  = fopen("g.dat", "w");
    fp_mat = fopen("mat.dat", "w");

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

	alphad = 5;           //AoA
	subalphad = 15;		  //subAoA
	leftright = 1.5;        //L
	updown = 2;           //H

    pi               = acos(-1.0); //PI

    alpha            = alphad * pi / 180.0; //
	subalpha		 = subalphad * pi / 180.0;

    cosalf           = cos(alpha);
    sinalf           = sin(alpha);

	subcosalf           = cos(subalpha);
    subsinalf           = sin(subalpha);

	//..allocate memory..
	xv  = new double[2 * nt];   zv = new double[2 * nt];
	xc  = new double[2 * nt];   zc  = new double[2 * nt];

	subxv  = new double[nt];   subzv = new double[nt];
	subxc  = new double[nt];   subzc  = new double[nt]; 
	
	pxv  = new double[2 * nt];   pzv = new double[2 * nt];
	pxc  = new double[2 * nt];   pzc  = new double[2 * nt];

	psubxv  = new double[nt];   psubzv = new double[nt];
	psubxc  = new double[nt];   psubzc  = new double[nt]; 

	xb = new double[nt+1];	sxb = new double[nt+1]; 
	subxb = new double[nt+1];	subsxb = new double[nt+1];

	psq = new double[kam2]; cam = new double[kam2]; pcam = new double[kam2];
	subpsq = new double[kam2]; subcam= new double[kam2]; psubcam= new double[kam2];

	int mxspl = (kam2-1)*4;

	cctk= new double[mxspl];cccm= new double[mxspl];enx= new double[2 * nt];
	enz = new double[2 * nt];   sig=  new double[2*nt];   gam= new double[2*nt]; pgam= new double[2*nt];

	subcctk= new double[mxspl];subcccm= new double[mxspl];subenx= new double[nt];
	subenz = new double[nt];   subsig=  new double[nt];   subgam= new double[nt]; subgam= new double[nt];

    return;
}

//-------------------------------------------------------------------
//   geom(): compute the location of vortex/source & control pts.
//-------------------------------------------------------------------
void geom()
{
    int i;
	double pper[kam2];
	double psubper[kam2];
	double p=0.5, m=0;
    double x_tilda;
    double *zup = new double[kam2];
    double *zlow= new double[kam2];

	double *pzup = new double[kam2];
    double *pzlow= new double[kam2];

	double *subzup = new double[kam2];
    double *subzlow= new double[kam2];

	double *psubzup = new double[kam2];
    double *psubzlow= new double[kam2];

	double thi[kam2];
	double subthi[kam2];

	double rx = 1.0;
	double ry = 0.0;
	double subrx = rx + leftright;
	double subry = updown;

	double subper[kam2];

	for (i = 0; i < kam2; i++)
	{
		subper[i] = per[i] + leftright;
	}

	for (i=0; i<kam2; i++)
	{
		thi[i]=(tc/0.2)*(0.2969*pow(per[i],0.5)-0.126*per[i]-0.3537*per[i]*per[i]
			+0.2843*pow(per[i],3.)-0.1015*pow(per[i],4.));

		subthi[i] = thi[i];
	}
    for (i=0; i<kam2; i++) {
        // Parabolic meanline: y=4*x(1-x)
        if(per[i]<=p)
			cam[i] = m * ( 2 * p * per[i] - per[i] * per[i] ) / ( p * p );
		else
			cam[i] = m * ( ( 1 - 2 * p ) + 2 * p * per[i]  - per[i] * per[i] ) /  ( 1 - 2 * p +  p * p );
  //     printf("%8.6f\n", cam[i]);
		subcam[i] = cam[i] + leftright;
        // Leading edge part is expanded by SQRT transformation
        // for interpolation of thickness 
        psq[i]       = sqrt(per[i]);
		subpsq[i]    = psq[i] + leftright;
    }

    uglydk(kam2, 1, 1, psq, thi, 0.0, 0.0, cctk);
    uglydk(kam2, 1, 1, per, cam, 0.0, 0.0, cccm);

	uglydk(kam2, 1, 1, subpsq, subthi, 0.0, 0.0, subcctk);
    uglydk(kam2, 1, 1, subper, subcam, 0.0, 0.0, subcccm);
	
    if (ispace == 0) {  // uniform spacing
       for (i=0; i<nt; i++) {
           xv[i]        = (double(i) + 0.25) / double(nt);
           xb[i]        = double(i) / double(nt);
           sxb[i]       = sqrt(xb[i]);
           xc[i]        = xv[i] + 0.5 * ( 1.0 / double(nt) );

		   subxv[i] = xv[i] + leftright;
		   subxb[i] = xb[i] + leftright;
		   subsxb[i] = sxb[i] + leftright;
		   subxc[i] = xc[i] + leftright;
		  
		   
	    }
           xb[nt]           = 1.0;
           sxb[nt]          = 1.0;
		   subxb[nt]        = 1.0 + leftright;
		   subsxb[nt]       = 1.0 + leftright;
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
	evaldk(kam2, nt, subper, subxv, subzv, subcccm);

    for (i=0; i<nt; i++ ) {               // vortex
        zv[i]       = zv[i] * fc;
		subzv[i]    = subzv[i] * fc + updown;
    }

    evaldk(kam2, nt, per, xc, zc, cccm);  // control point
	evaldk(kam2, nt, subper, subxc, subzc, subcccm);

    for (i=0; i<nt; i++) {
        zc[i]       = zc[i] * fc;
		subzc[i]    = subzc[i] * fc + updown;
    }

	////////////////////////////    mian   ///////////////////////////
	fprintf(fp_geo, "variables = x, z\n");
    fprintf(fp_geo, "zone i = %d\n", nt);
	for (i=0; i<nt; i++){	
		pxv[i] = ( xv[i] - rx ) * realcos( -alphad ) - ( zv[i] - ry ) * realsin( -alphad ) + rx;
		pzv[i] = ( xv[i] - rx ) * realsin( -alphad ) - ( zv[i] - ry ) * realcos( -alphad ) + ry;
		fprintf(fp_geo, "%0.5lf  %0.5lf\n",pxv[i],pzv[i]);
	}
	fprintf(fp_geo, "zone i = %d\n", nt);
	for (i=0; i<nt; i++){
		pxc[i] = ( xc[i] - rx ) * realcos( -alphad ) - ( zc[i] - ry ) * realsin( -alphad ) + rx;
		pzc[i] = ( xc[i] - rx ) * realsin( -alphad ) - ( zc[i] - ry  ) * realcos( -alphad ) + ry;
		fprintf(fp_geo, "%0.5lf  %0.5lf\n",pxc[i],pzc[i]);
	}

    for (i=0; i<kam2; i++) {
        zup[i]       = cam[i] * fc + thi[i];
        zlow[i]      = cam[i] * fc - thi[i];
    }
	
	for (i=0; i<kam2; i++){
		pper[i] = ( per[i] - rx ) * realcos( -alphad ) - ( zup[i] - ry ) * realsin( -alphad ) + rx;
		pzup[i] = ( per[i] - rx ) * realsin( -alphad ) - ( zup[i] - ry ) * realcos( -alphad ) + ry;
		pzlow[i] = ( per[i] - rx ) * realsin( -alphad ) - ( zlow[i] - ry ) * realcos( -alphad ) + ry;
		pcam[i] = ( per[i] - rx ) * realsin( -alphad ) - ( cam[i] - ry ) * realcos( -alphad ) + ry;
	}

    fprintf(fp_geo, "variables = x, z\n");
    fprintf(fp_geo, "zone i = %d\n", kam2);
    for (i=0; i<kam2; i++) {
        fprintf(fp_geo,"%f %f\n", pper[i], pzup[i]);
    }
    fprintf(fp_geo, "zone i = %d\n", kam2);
    for (i=0; i<kam2; i++) {
        fprintf(fp_geo,"%f %f\n", pper[i], pzlow[i]);
    }
//    fprintf(fp_geo, "zone i = %d\n", kam2);
//    for (i=0; i<kam2; i++) {
//        fprintf(fp_geo,"%f %f\n", pper[i], pcam[i]*fc);
//    }
	////////////////////////////    mian   ///////////////////////////

	////////////////////////////    sub   ///////////////////////////
	fprintf(fp_geo, "variables = x, z\n");
    fprintf(fp_geo, "zone i = %d\n", nt);
	for (i=0; i<nt; i++){
		psubxv[i] = ( subxv[i] - subrx ) * realcos( -subalphad ) - ( subzv[i] - subry ) * realsin( -subalphad ) + subrx;
		psubzv[i] = ( subxv[i] - subrx ) * realsin( -subalphad ) - ( subzv[i] - subry ) * realcos( -subalphad ) + subry;
		fprintf(fp_geo, "%0.5lf  %0.5lf\n",psubxv[i],psubzv[i]);
	}

	fprintf(fp_geo, "zone i = %d\n", nt);
	for (i=0; i<nt; i++){
		psubxc[i] = ( subxc[i] - subrx ) * realcos( -subalphad ) - ( subzc[i] - subry ) * realsin( -subalphad ) + subrx;
		psubzc[i] = ( subxc[i] - subrx ) * realsin( -subalphad ) - ( subzc[i] - subry ) * realcos( -subalphad ) + subry;

		fprintf(fp_geo, "%0.5lf  %0.5lf\n",psubxc[i],psubzc[i]);
	}
	
    for (i=0; i<kam2; i++) {
        subzup[i]       = subcam[i] * fc + subthi[i] + updown;
        subzlow[i]      = subcam[i] * fc - subthi[i] + updown;
    }

	for (i=0; i<kam2; i++){
		psubper[i] = ( subper[i] - subrx ) * realcos( -subalphad ) - ( subzup[i] - subry ) * realsin( -subalphad ) + subrx;
		psubzup[i] = ( subper[i] - subrx ) * realsin( -subalphad ) - ( subzup[i] - subry ) * realcos( -subalphad ) + subry;
		psubzlow[i] = ( subper[i] - subrx ) * realsin( -subalphad ) - ( subzlow[i] - subry ) * realcos( -subalphad ) + subry;
		psubcam[i] = ( subper[i] - subrx ) * realsin( -subalphad ) - ( subcam[i] - subry ) * realcos( -subalphad ) + subry;
	}

    fprintf(fp_geo, "variables = x, z\n");
    fprintf(fp_geo, "zone i = %d\n", kam2);
    for (i=0; i<kam2; i++) {
        fprintf(fp_geo,"%f %f\n", psubper[i], psubzup[i]);
    }
    fprintf(fp_geo, "zone i = %d\n", kam2);
    for (i=0; i<kam2; i++) {
        fprintf(fp_geo,"%f %f\n", psubper[i], psubzlow[i]);
    }
//    fprintf(fp_geo, "zone i = %d\n", kam2);
//    for (i=0; i<kam2; i++) {
//        fprintf(fp_geo,"%f %f\n", psubper[i], psubcam[i]*fc);
//    }
	////////////////////////////    sub   ///////////////////////////
	
	for(i=0;i<nt;i++)
	{
		xv[i+nt]=subxv[i];
		zv[i+nt]=subzv[i];
		xc[i+nt]=subxc[i];
		zc[i+nt]=subzc[i];
		xb[i+nt]=subxb[i];
		sxb[i+nt]=subsxb[i];
	}


	for(i=0;i<kam2;i++)
	{
		per[i+kam2]=subper[i];
		psq[i+kam2]=subpsq[i];		
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


    double *tcn = new double[2 * (nt+1)];
    double *dc  = new double[2 * nt];
    double *d2c = new double[2 * nt];
    double *zvi = new double[2 * nt];
    double *ut  = new double[2 * nt];
    double *wt  = new double[2 * nt];

    //-- alloc memory for coefficient matrices ---
    ipvt= new int[2 * nt];  

    a    = new double*[2 * nt];
    for (i=0; i<2 * nt; i++) {
        a[i] = new double[2 * nt];
    }

    alud = new double*[2 * nt];
    for (i=0; i<2 * nt; i++) {
        alud[i] = new double[2 * nt];
    }

    if ( ispace == 0 ) {
        evaldk(2 * kam2, 2*(nt+1), psq, sxb, tcn, cctk);
    }
    else {
        evaldk(2 * kam2, 2*(nt+1), per, xb, tcn, cctk);
    }

    for (i=0; i<2 * nt; i++) {
        sig[i]       = tc * (tcn[i+1] - tcn[i]);
    }

    drivdk(2 * kam2, 2 * nt, per, xc, dc, d2c, cccm);

    for (i=0; i<2 * nt; i++) {
        slope        = fc * dc[i];
        theta        = atan(slope);
        enx[i]       = -sin(theta);
        enz[i]       = cos(theta);
    }

    dtp              = 0.5 / pi;
    for (i=0; i<2 * nt; i++) {
        ut[i]        = 0.0;
        wt[i]        = 0.0;
        for (j=0; j<2 * nt; j++) {
            dx       = xc[i] - xv[j]; //delta x
            dz       = zc[i] - zv[j]; //delta z
            dtpr    = dtp / (dx*dx + dz*dz); //r
            unitv    = dtpr * dz; //uxi
            wnitv    = -dtpr * dx; //wij
            units    = dtpr * dx;
            wnits    = dtpr * dz;
            a[i][j]  = enx[i] * unitv + enz[i] * wnitv; //enx:nxi  enz:nzi
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
	for (i=nt; i<2 * nt; i++) {
        gam[i] = -enx[i] * (subcosalf+ut[i]) - enz[i] * (subsinalf+wt[i]);
        fprintf(fp_mat, "%f   \n", gam[i]);
    }

    //--- LU decomposition of coeff. matrix ---
    ludcmp(a, alud, 2 * nt, ipvt, &ludsign);
    lubksb(alud, 2 * nt, ipvt, gam);

    fprintf(fp_gam, "variables = x, y\n");
    fprintf(fp_gam, "zone i = %d\n",nt);
    for (i=0; i<nt; i++) {
        fprintf(fp_gam, "%f %f\n", xc[i], gam[i]);
    }
	fprintf(fp_gam, "zone i = %d\n",nt);
    for (i=nt; i<2 * nt; i++) {
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
    double tgam, dx, dz, bug, fxtt, fztt, fxlt, fzlt, force, cl, subforce, subcl, subfxtt,
		subfztt, subfxlt, subfzlt, velup, velow, delnt, subvelup, subvelow, subdelnt;
    double *voncom=new double[2 * nt];   double *vind = new double[2 * nt];
    double *cpup  =new double[2 * nt];   double *cplow= new double[2 * nt];
    double *cpnega=new double[2 * nt];   double *dgam = new double[2 * nt];
    double *uts  = new double[2 * nt];   double *wts  = new double[2 * nt];
    double *uls  = new double[2 * nt];   double *wls  = new double[2 * nt];
    double *fxt  = new double[2 * nt];   double *fzt  = new double[2 * nt];
    double *fxl  = new double[2 * nt];   double *fzl  = new double[2 * nt];

    tgam             = 0.0;
    for (i=0; i<2 * nt; i++) {
        tgam        += gam[i];
        uts[i]       = 0.0;
        wts[i]       = 0.0;
        uls[i]       = 0.0;
        wls[i]       = 0.0;
        for (j=0; j<2 * nt; j++) {
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

	subfxtt             = 0.0;
    subfztt             = 0.0;
    subfxlt             = 0.0;
    subfzlt             = 0.0;

    for (i=0; i< nt; i++) {
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

	for (i=nt; i< 2 * nt; i++) {
        fxt[i]       = -sig[i] * (subcosalf + uts[i] + uls[i]);
        fzt[i]       = -sig[i] * (subsinalf + wts[i] + wls[i]);
        fxl[i]       = -gam[i] * (subsinalf + wts[i] + wls[i]);
        fzl[i]       = gam[i] * (subcosalf + uts[i] + uls[i]);
        subfxtt        += fxt[i];
        subfztt        += fzt[i];
        subfxlt        += fxl[i];
        subfzlt        += fzl[i];
    }
    subforce            = sqrt(pow((subfxtt+subfxlt),2) + pow((subfztt+subfzlt),2));
    subcl               = 2.0 * subforce;

	displayright();
    printf("Lift Ceofficient (2-D VLM) of first  =  %f\n", cl);
	printf("Lift Ceofficient (2-D VLM) of second  =  %f\n", subcl);
	printf("Lift Ceofficient(2*pi*sin(alpha)) =  %f\n", 2*pi*sin(alpha));
	printf("Lift Ceofficient(4*pi*f/c) =  %f\n", 4*pi*fc);


    delnt            = 1.0 / double(nt);
	subdelnt         = 1.0 / double(nt);
    for (i=0; i< nt; i++) {
        voncom[i]    = cosalf * enz[i] - sinalf * enx[i];
        vind[i]      = gam[i] / (2.0 * double(delnt));
        velup        = voncom[i] + vind[i];
        velow        = voncom[i] - vind[i];
        cpup[i]      = 1.0 -  pow(velup,2);
        cplow[i]     = 1.0 -  pow(velow,2);
        cpnega[i]    = cpup[i] - cplow[i];
        dgam[i]      = gam[i] / double(delnt);
    }
	for (i=nt; i< 2 * nt; i++) {
        voncom[i]    = subcosalf * enz[i] - subsinalf * enx[i];
        vind[i]      = gam[i] / (2.0 * double(delnt));
        subvelup        = voncom[i] + vind[i];
        subvelow        = voncom[i] - vind[i];
        cpup[i]      = 1.0 -  pow(subvelup,2);
        cplow[i]     = 1.0 -  pow(subvelow,2);
        cpnega[i]    = cpup[i] - cplow[i];
        dgam[i]      = gam[i] / double(subdelnt);
    }
    
	printf("Lift Ceofficient (2*Gam/Uc)   =  %f\n", 2*tgam);
    fprintf(fp_pre, "variables = x, z\n");
    fprintf(fp_pre, "zone i = %d\n", nt);
    for (i=0; i<nt; i++) {
        fprintf(fp_pre,"%f %f\n", xc[i], -cpup[i]);
    }

	fprintf(fp_pre, "zone i = %d\n", nt);
	for (i=nt; i<2 * nt; i++) {
       fprintf(fp_pre,"%f %f\n", xc[i], -cpup[i]);
    }

    fprintf(fp_pre, "zone i = %d\n", nt);
    for (i=0; i<nt; i++) {
        fprintf(fp_pre,"%f %f\n", xc[i], -cplow[i]);
    }
	fprintf(fp_pre, "zone i = %d\n", nt);
    for (i=nt; i<2 * nt; i++) {
        fprintf(fp_pre,"%f %f\n", xc[i], -cplow[i]);
    }
    return;
}

double realsin( double angle )
{
	double radian;
	radian = angle * PI / 180;
	return sin( radian );
}
double realcos( double angle )
{
	double radian;
	radian = angle * PI / 180;
	return cos( radian );
}
