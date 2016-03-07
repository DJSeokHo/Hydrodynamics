/////////////////////////////////////////////////////////////////
// PropDLL: Propeller Dyanmic Link Library
//
/////////////////////////////////////////////////////////////////
// Contact: Prof. Lee, Chang-Sup 
//          Dept of Naval Architecture and Ocean Engineering
//          Chungnam National University
//          Voice: 82-42-821-6623
//          Fax: 82-42-823-8763
//          E-mail: cslee@naoe.chungnam.ac.kr
//
/////////////////////////////////////////////////////////////////
// PropDLL contains duck routines, tension spline routines and 
//         linear equation solvers.
//
// References:
//   (1) Ducklib Functions:                             
//       Kerwin, J.E., M.I.T.                           
//   (2) Tension Spline Functions:                      
//       McCartin, B. J., "Application of exponential   
//          spline in computational fluid dynamics,"    
//          AIAA J., 1983                               
//   (3) Linsol Functions:                              
//       Press et al, "Numerical Recipes in C", 1988    
//================================================================
// Release Version 1.0, April 15, 1996, LCS
//
// Language: C/C++
//
// Usage: (1) Insert #include "c:\lib\proplib.h" in the source file
//        (2) Include the directory c:\lib\PropDLL.lib in Link 
//            option
//        (3) PropDLL.DLL should reside in c:\windows\system
//        (4) PropDLL.lib should reside in c:\lib
//        (5) proplib.h should reside in c:\lib
/////////////////////////////////////////////////////////////////
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <malloc.h>
#include <math.h>
#include "proplib.h"                  

using namespace std;

ofstream fo_err("proplib.err", ios::in);
///////////////////////////////////////////////////////////////////

//********************************************************************
void linsol(double **a, double **alud, int n, int *indx, double *b,
            double *x, int kase, int inital, int itmax, double tolx)
{
int i;
double d, tiny;

/*       This routine calls some useful routines given by the text
 *     "NUMERICAL RECIPES in C", to solve the set of N linear equations
 *     A.X=B by LU decomposition either with or without iterative
 *     improvement of a solution.
 *       A, B, N, KASE and INITAL are input arguments.
 *       ALUD and INDX are either output or input, depending on KASE.
 *     Here ALUD[0,...,n-1][0,...,n-1] is, the LU decomposition of the
 *     matrix A[0,...,n-1][0,...,n-1], determined by the routine LUDCMP.
 *     INDX[0,...,n-1] is a vector which records the row permutation
 *     effected by the partial pivoting.
 *       X[0,...,n-1] is finally an output solution vector and also if
 *     KASE=3 or 4 with INITAL=1, it is input as an initial guess of
 *     solution supplied by a user. In this case, ITMAX (maximum no. of
 *     iterations allowed) and TOLX (tolerance used in convergence test)
 *     are selected by a user.
 *       If KASE=1 or 3, this routine calls the routine LUDCMP to provide
 *     an output matrix ALUD and an output vector INDX  determined by
 *     the LU decomposition of A, otherwise ALUD and INDX should be
 *     entered as input arguments returned already from LUDCMP.
 *       INITAL is a control data of either 1 or 0, depending on whether
 *     the initial guess for solution is given by user or not.
 *
 *       If you want to solves the linear set of equations by performing
 *     LU decomposition directly, you enter KASE=1;
 *       If you subsequently want to solve a set of equations with the
 *     same matrix A but a different right-hand side B, you enter KASE=2
 *     so that only the routine LUBKSB is done. Of course the matrix ALUD
 *     and INDX were already returned from LUDCMP;
 *       If you want to solve large sets of linear equations, it is
 *     recommended to enter KASE=3 in order to use iterated improvement
 *     of the solution with LU composition procedure.
 *       If you already have the LU decomposition form of the matrix A,
 *     you enter KASE=4.
 *
 *       A, ALUD, N, INDX and B  are not modified by this routine and
 *     can be left in place for successive calls with different right-hand
 *     sides B.
 *
 *     Note: 1) C_version 1.0, translated from KPA4a, Jan. 7, 1992, LCS
 *           2) Note differences in arguments with KPA4a.
 *           3) C_version 1.1, zero-offset implemented, Mar. 17, 1993, LCS
 */
		
    if ((kase<1) || (kase>4)) {
        fo_err << "Error: not assigned cases in LINSOL" << endl;
        fo_err.close();
        // abort; **** kase: should be checked before calling LINSOL
    }

    /* direct LU decomposition method */
    if (kase==1) {
        ludcmp(a, alud, n, indx, &d);
        for (i=0; i<n; i++) {
                x[i] = b[i];
        }
        lubksb(alud, n, indx, x);
    }

    /* subsequent use of LU-method with the same A and a different B */
    if (kase==2) {
        for (i=0; i<n; i++) {
            x[i] = b[i];
        }
        lubksb(alud, n, indx, x);
    }

    /* iterated improvement with LU decomposition */
    if (kase==3) {
        ludcmp(a, alud, n, indx, &d);
        tiny = 1.0e-20;
        if (inital==0) {
            for (i=0; i<n; i++) {
                if (fabs(a[i][i])<tiny) {
                    x[i] = b[i]/tiny;
                }
                else {
                    x[i] = b[i]/a[i][i];
                }
            }
        }
       iterlu(a, alud, n, indx, b, x, itmax, tolx);
    }

    /* subsequent use of iterated improvement with same A and a diff. B */
    if (kase==4) {
        tiny = 1.0e-20;
        if (inital==0) {
            for (i=0; i<n; i++) {
                if( fabs(a[i][i])<tiny ) {
                    x[i] = b[i]/tiny;
                }
                else {
                    x[i] = b[i]/a[i][i];
                }
            }
        }
        iterlu(a, alud, n, indx, b, x, itmax, tolx);
    }

    return;
}
/**********************************************************************/

#define TINY 1.0e-20;
void ludcmp(double **a,double **alud,int n,int *indx,double *d)
{
        int i,imax,j,k;
        double *vv,big,dum,sum,temp;
/*      double *vector(),**matrix();
        void nrerror(),free_vector();        */

/*      The routine given on page 43, 'NUMERICAL RECIPES in C', is slightly
 *    modified to take our purpose. Given an N x N matrix A,
 *    A[0,...,n-1][0,...,n-1],
 *    this routine provides ALUD by the LU decomposition of a rowwise
 *    permutation  of itself.
 *      A and N are input, ALUD is output, arranged as in equation
 *    (2.3.14) of 'NUMERICAL RECIPES'.
 *    INDX[0,...,n-1] is an output vector which
 *    records the row permutation effected by the partial pivoting.
 *    d is output as +1 or -1 depending on number of interchange.
 *      This routine is used in combination with LUBKSB to solve linear
 *    equations.
 */
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++)
            alud[i][j]=a[i][j];
    }

    vv=vector(0,n-1);
    *d=1.0;
    for (i=0; i<n; i++) {
        big=0.0;
        for (j=0; j<n; j++)
            if ((temp=fabs(alud[i][j])) > big) big=temp;
        if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
        vv[i]=1.0/big;
    }
    for (j=0; j<n; j++) {
        for (i=0; i<j; i++) {
            sum=alud[i][j];
            for (k=0; k<i; k++) sum -= alud[i][k]*alud[k][j];
            alud[i][j]=sum;
        }
        big=0.0;
        for (i=j; i<n; i++) {
            sum=alud[i][j];
            for (k=0; k<j; k++)
                sum -= alud[i][k]*alud[k][j];
            alud[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=0; k<n; k++) {
                dum=alud[imax][k];
                alud[imax][k]=alud[j][k];
                alud[j][k]=dum;
            }
            *d = -(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (alud[j][j] == 0.0) alud[j][j]=TINY;
        if (j != n-1) {
            dum=1.0/(alud[j][j]);
            for (i=j+1; i<n; i++) alud[i][j] *= dum;
        }
    }
    free_vector(vv,0,n-1);
}

#undef TINY
/**********************************************************************/
void lubksb(double **alud, int n, int *indx, double *b)
{
    int i,ii=0,ip,j,non_vanishing_index;
    double sum;

/*       This routine is taken from pages 44, 'NUMERICAL RECIPES in C',
 *     with change of variable A into ALUD. This routine solves the set
 *     of N linear equations A.X=B. Here ALUD is input with
 *     ALUD[0,...,n-1][0,...,n-1],
 *     the LU decomposition of original A matrix, determined by the routine
 *     LUDCMP. INDX[0,...,n-1] is input as the permutation vector returned by
 *     LUDCMP. B[1..n] is input as the right-hand side vector and return
 *     with the solution vector X[0,...,n-1].
 *       ALUD, N and INDX are not modified by this routine and can be
 *     left in place for successive calls with different right-hand sides B.
 */
    for (i=0; i<n; i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii) {
            for (j=non_vanishing_index; j<=i-1; j++)
                sum -= alud[i][j]*b[j];
        } else if (sum) {
            ii=1; non_vanishing_index = i;
        }
        b[i]=sum;
    }
    for (i=n-1; i>=0; i--) {
        sum=b[i];
        for (j=i+1; j<n; j++) sum -= alud[i][j]*b[j];
        b[i]=sum/alud[i][i];
    }
}
/**********************************************************************/
void iterlu(double **a, double **alud, int n, int *indx, double *b,
            double *x, int itmax, double tolx)
{
    int j,i,it;
    double sdp;
    double *r,rmax;
/*  double *vector();
    void lubksb(),free_vector();       */

/*     This routine is taken from page 50, 'NUMERICAL RECIPES in C', with
 *     adding several times of iterations. This routine solves the set of
 *     N linear equations A.X=B by iterative mothod.
 *       The matrix A[0,...,n-1][0,...,n-1], and vectors B[0,...,n-1] and
 *     X[0,...,n-1] are input.      The matrix ALUD[0,...,n-1][0,...,n-1]
 *     is also input as its LU decomposition of A, returned by the routine
 *     LUDCMP. INDX[0,...,n-1] is input as the permutation vector returned
 *     by LUDCMP. B[0,...,n-1] is input as the right-hand side vector.
 *       Only the solution vector X[0,...,n-1] is modified to an improved
 *     set of values. ITMAX is input as no. of maximum iterations allowed.
 *     TOLX is input as tolerance of used in convergence test.
 */
    r=vector(0,n-1);

    for (it=0; it<itmax; it++){
        for (i=0; i<n; i++) {
            sdp = -b[i];
            for (j=0; j<n; j++) sdp += a[i][j]*x[j];
            r[i]=sdp;
        }
        lubksb(alud, n, indx, r);
        for (i=0; i<n; i++) x[i] -= r[i];
        rmax=0.0;
        for (i=0; i<n; i++){
            if (fabs(r[i]) > rmax)
                rmax = fabs(r[i]);
        }
        if (rmax < tolx) {
            fo_err << " ITERLU: No. of iteration = " << it << endl;
            return;
        }
    }
    fo_err << "No convergence of iterative method of solution" << endl;
    fo_err << " ITERLU: No. of iteration = " << it-1 << endl;
    fo_err << " ITERLU: Max. residual    = " << rmax << endl;
    fo_err.close();
    free_vector(r,0,n-1);
    return; /// abort; **** Add iteration number check routine
}
    
/**********************************************************************/
/*--------------------Numerical Recipes Routines----------------------*/
void nrerror(char error_text[])
{

    fo_err <<"Numerical Recipes run-time error..." << endl;
    fo_err << error_text << endl ;
    fo_err <<"...now exiting to system..." << endl;
	  fo_err.close();
    abort();
}

double *vector(int nl, int nh)
{
    double *v;

    v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl;
}

int *ivector(int nl, int nh)
{
    int *v;

    v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
    if (!v) nrerror("allocation failure in ivector()");
    return v-nl;
}

double **matrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;

    m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m -= nrl;

    for(i=nrl; i<=nrh; i++) {
        m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
        if (!m[i]) nrerror("allocation failure 2 in matrix()");
        m[i] -= ncl;
    }
    return m;
}


void free_vector(double *v, int nl, int nh)
{
    free((char*) (v+nl));
}

void free_ivector(int *v, int nl, int nh)
{
    free((char*) (v+nl));
}

void free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
    int i;

    for(i=nrh; i>=nrl; i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}

/**********************************************************************/
/*-------------------- SPLINE.C ---------------------------------------
*
*     SPLINE.C
*             Modified from 1975 DUCK series  by  J.E.Kerwin
*             C-Version 2.1, Jan. 7, 1992, LCS
*                       2.2, Mar. 17, 1993, LCS
*             Note: 1) Zero-offset convention is implemented in all
*                      duck library routines.
* ----------------------------------------------------------------------
*/

/*---------------------------------------------------------------------
*    uglydk: prepares the coeff. matrix for cubic spline.
*----------------------------------------------------------------------
*/
void uglydk(int nin,int ncl,int ncr, double *Xin, double *Yin,
            double esl,double esr, double *Ae)
{
    //bug(kgd)2003.6.26..double h[50],d[50],a[2500],s[50],slp,hold,bug;
    double h[1000],d[1000],a[100000],s[1000],slp,hold,bug;
    double  half=0.5, two=2.0, six=6.0, rad=1.745329e-02;
    int     n,nm1,nm2,nm3, j,l,m,neq,nsq,kerror;

    /*....Zero-offset conversion...3/17/1993...*/
    double *xin, *yin, *ae, *a_simq, *s_simq;
    xin = Xin - 1;
    yin = Yin - 1;
    ae  = Ae - 1;
    a_simq = a + 1;
    s_simq = s + 1;

    nm1=nin-1;
    nm2=nm1-1;
    nm3=nm2-1;
    neq=nm2;
    for (n=1; n<=nm1; n++) {
        h[n] = xin[n+1] -xin[n];
        d[n] = ( yin[n+1] - yin[n] ) / h[n];
    }

    if (ncl == 2) neq += 1;
    if (ncr == 2) neq += 1;
    nsq = neq * neq;

    for (n=1; n<=nsq; n++) a[n] = 0.0;

    j=1;   l=1;

    if (!(ncl<2)) {
        a[1] = two * h[1];
        a[2] = h[1];
        slp  = esl * rad;
        s[1] = ( d[1] - tan(slp) ) * six;
        j += 1;
        l += neq+1;
        a[l-1] = h[1];
    }

    for (n=1; n<=nm2; n++,j++,l+=neq+1) {
        if (n > 1)  a[l-1] = h[n];
        a[l] = two * ( h[n] + h[n+1] );
        if (n < nm2) a[l+1] = h[n+1];
        if (n == 2 && ncl == 1) a[l-1] -= h[n]*h[n]/h[n+1];
        if (n == 1 && ncl == 1) a[l]   += (1.0+h[n]/h[n+1])*h[n];
        if (n == nm2 && ncr == 1) a[l] += (1.0+h[n+1]/h[n])*h[n+1];
        if (n == nm3 && ncr == 1) a[l+1]-= h[n+2]*h[n+2]/h[n+1];
        s[j] = (d[n+1] - d[n]) * six;
    }

    if (!(ncr<2))  {
        a[l-1] = h[nm1];
        a[l]   = -two*h[nm1];
        l     -= neq;
        a[l]   = -h[nm1];
        slp    = esr*rad;
        s[j]   = (d[nm1] + tan(slp))*six;
    }

    simq(a_simq,s_simq,neq,&kerror);

    hold = s[neq];

    if (!(ncl==2))  {
        for (n=1; n<=nm2; n++) {
            m    = nm2-n+2;
            s[m] = s[m-1];
        }
        if (ncl == 0) s[1] = 0.0;
        bug = h[1] / h[2];
        if (ncl == 1) s[1] = (1.0+bug)*s[2]-bug*s[3];
    }

    if (ncr == 0) s[nin] = 0.0;

    bug = h[nm1] / h[nm2];
    if (ncr == 1) s[nin] = (1.0+bug)*s[nm1]-bug*s[nm2];
    if (ncr == 2) s[nin] = hold;

    for (n=1; n<=nm1; n++) {
        ae[n] = (s[n+1]-s[n]) / (six*h[n]);
        m     = n + nm1;
        ae[m] = half * s[n];
        m    += nm1;
        ae[m] = d[n] - h[n] * (two*s[n]+s[n+1]) / six;
        m    += nm1;
        ae[m] = yin[n];
    }
}

/*---------------------------------------------------------------------
*    evaldk: interpolates by splines. uses coeff. matrix made by uglydk.
*----------------------------------------------------------------------
*/
void evaldk(int nin, int nout,   double *Xin, double *Xout,
            double *Yout, double *A)
{
    int j,n,nm1,mout,j2,j3,j4;//..j1
    double del,h1,h2,h3;

    /*....Zero-offset conversion...3/17/1993...*/
    double *xin, *xout, *yout, *a;
    xin  = Xin -1;
    xout = Xout -1;
    yout = Yout -1;
    a    = A -1;

    nm1  = nin - 1;
    mout = abs(nout);

    if (nout <= 0) {
        del  = (xin[nin] - xin[1]) / (mout-1);
        for (n=1; n<=mout; n++) {
            xout[n]  = xin[1] + (n-1)*del;
        }
    }

    j  = 1;
    for (n=1; n<=mout; n++) {
        if (xout[n] >= xin[2]  &&  xout[n] < xin[nm1]) {
            while ( xout[n] >= xin[j+1] ) j++;
            while ( xout[n] <  xin[j]   ) j--;
        }
        else {
            if( xout[n] <  xin[2]   ) j=1;
            if( xout[n] >= xin[nm1] ) j=nm1;
        }

        h1  = xout[n] - xin[j];
        h2  = h1 * h1;
        h3  = h1 * h2;
        j2  = j  + nm1;
        j3  = j2 + nm1;
        j4  = j3 + nm1;
        yout[n] = a[j]*h3 + a[j2]*h2 + a[j3]*h1 + a[j4];
    }
}


/*---------------------------------------------------------------------
*    drivdk: computes the first and second derivatives by spline.
*----------------------------------------------------------------------
*/
void drivdk(int nin, int nout, double *Xin, double *Xout,
            double *Dydx, double *D2ydx, double *A)
{
    int j,n,nm1,j2,j3;//..j1
    double h1,h2;

    /*....Zero-offset conversion...3/17/1993...*/
    double *xin, *xout, *dydx, *d2ydx, *a;
    xin   = Xin -1;
    xout  = Xout -1;
    dydx  = Dydx -1;
    d2ydx = D2ydx - 1;
    a     = A - 1;

    nm1  = nin - 1;
    j    = 1;

    for (n=1; n<=nout; n++) {
        if (xout[n] >= xin[2]  &&  xout[n] < xin[nm1]) {
            while (xout[n] >= xin[j+1]) j++;
        }
        else {
            if (xout[n] <  xin[2]  ) j=1;
            if (xout[n] >= xin[nm1]) j=nm1;
        }

        h1       = xout[n] - xin[j];
        h2       = h1 * h1;
        j2       = j  + nm1;
        j3       = j2 + nm1;
        dydx[n]  = 3.0*a[j]*h2 + 2.0*a[j2]*h1 + a[j3];
        d2ydx[n] = 6.0*a[j]*h1 + 2.0*a[j2];
    }
}

/*---------------------------------------------------------------------
*    intedk: evaluates area and the 1st/2nd moments of area by spline.
*----------------------------------------------------------------------
*/
/*
void intedk(int nin, double *Xin, double xl, double xu,
            double &ydx, double &xydx, double &xxydx, double *A)
*/
void intedk(int nin, double *Xin, double xl, double xu,
            double *ydx, double *xydx, double *xxydx, double *A)
{
    int n,jl,ju,j1,j2,j3,j4,nm1;
    double h1,h2,h3,h4,h5,h6,bug,cat,dog,pig;

    /*....Zero-offset conversion...3/17/1993...*/
    double *xin, *a;
    xin   = Xin - 1;
    a     = A - 1;

    nm1   = nin - 1;
    /* determine lower bound index  jl */
    if (xl <  xin[1]) {
        jl = 1;
    }
    else {
        jl = nm1;
        for (n=1; n<=nin; n++) {
            if (xl < xin[n]) {
                jl  = n - 1; break;
            }
        }
    }
    /* determine upper bound index  ju */
    if (xu >= xin[nin]) {
        ju = nm1;
    }
    else {
		if (xu <= xin[1]) {
			ju = 1;
		} else {
           ju = 1;
           for (n=jl; n<=nin; n++) {
               if (xu < xin[n]) {
                  ju  = n - 1; break;
               }
           }
		}
    }
/* print out  matrix coeff.
    printf(" lower index = %d, upper index = %d\n", jl,ju);
    for (n=1; n<=nin*nin-4; n++) {
        printf(" ae[%d]= %9f \n", n,a[n]);
    }
*/
    h1  = xl -xin[jl];
    h2  = h1 * h1;
    h3  = h1 * h2;
    h4  = h2 * h2;
    h5  = h2 * h3;
    h6  = h3 * h3;
    j1  = jl;
    j2  = jl + nm1;
    j3  = j2 + nm1;
    j4  = j3 + nm1;

    *ydx   = -a[j1]/4.0*h4 - a[j2]/3.0*h3 - a[j3]/2.0*h2 - a[j4]*h1;
    bug    = -a[j1]/5.0*h5 - a[j2]/4.0*h4 - a[j3]/3.0*h3 - a[j4]/2.0*h2;
    *xydx   = bug + xin[j1]*(*ydx);
    bug   = -a[j1]/6.0*h6 - a[j2]/5.0*h5 - a[j3]/4.0*h4 - a[j4]/3.0*h3;
    *xxydx = bug + 2.0*xin[j1]*(*xydx) - xin[j1]*xin[j1]*(*ydx);

    for (n=jl; n<=ju; n++) {
        h1    = xin[n+1] - xin[n];
        if( n == ju )   h1  = xu - xin[n];
        h2     = h1 * h1;
        h3     = h1 * h2;
        h4     = h2 * h2;
        h5     = h2 * h3;
        h6     = h3 * h3;
        j1     = n;
        j2     = j1 + nm1;
        j3     = j2 + nm1;
        j4     = j3 + nm1;
        bug    = a[j1]/4.0*h4 + a[j2]/3.0*h3 + a[j3]/2.0*h2 + a[j4]*h1;
        *ydx  += bug;
        cat    = a[j1]/5.0*h5 + a[j2]/4.0*h4 + a[j3]/3.0*h3 + a[j4]/2.0*h2;
        pig    = cat + xin[j1]*bug;
        *xydx += pig;
        dog    = a[j1]/6.0*h6 + a[j2]/5.0*h5 + a[j3]/4.0*h4 + a[j4]/3.0*h3;
        *xxydx+= dog + 2.0*xin[j1]*pig - xin[j1]*xin[j1]*bug;
    }
/*  printf(" Area = %6.3f, 1st Mt = %6.3f, 2nd Mt = %6.3f\n",
             *ydx,*xydx,*xxydx);  */
}


/*---------------------------------------------------------------------
*    simq: solves a linear system of equations.
*----------------------------------------------------------------------
*/
void simq(double *A, double *B, int n, int *ks)
{
    int j,jj,jjx,jx,jy,it,i,ij,ix,ixj,ixjx,imax,i1,i2,k,iqs;//..ixy
    int ny,ia,ib,ic;
    double tol=0.0,biga,save;

    /*....Zero-offset conversion...3/17/1993...*/
    double *a, *b;
    a = A - 1;
    b = B - 1;

/* Print the contents of matrix and the right-side vector */
/*  printf(" Content of coeff. matrix\n ");
    for (j=1; j<=n*n; j++) {
        printf(" j=%3d, a[%3d]=%9f\n",j,j,a[j]);
    }
    printf(" Right-hand side vector\n ");
    for (j=1; j<=n; j++) {
        printf(" j=%3d, b[%3d]=%9f\n",j,j,b[j]);
    }
*/

    *ks = 0;
    jj = -n;
/* ------------- beginning of big for-loop ------------- */
    for (j=1; j<=n; j++) {
        jy  = j + 1;
        jj += n + 1;
        biga= 0.0;
        it  = jj - j;
        for (i=j; i<=n; i++) {
            ij  = it + i;
            if (fabs(biga) < fabs(a[ij]) ) {
                biga = a[ij];
                imax = i;
            }
        }

        if (fabs(biga) <= tol) {  *ks = 1; return;  };

        i1 = j + n*(j-2);
        it = imax - j;
        for (k=j; k<=n; k++) {
            i1   += n;
            i2    = i1 + it;
            save  = a[i1];
            a[i1] = a[i2];
            a[i2] = save;
            a[i1] = a[i1] / biga;
        }

        save   = b[imax];
        b[imax]= b[j];
        b[j]   = save / biga;

        if (j == n) break;   /* break the for loop if j==n */

        iqs = n * (j-1);

        for (ix=jy; ix<=n; ix++) {
            ixj  =  iqs + ix;
            it   =  j   - ix;
            for (jx=jy; jx<=n; jx++) {
                ixjx  = n*(jx-1) + ix;
                jjx   = ixjx     + it;
                a[ixjx] -= a[ixj] * a[jjx];
            }
            b[ix]  -= b[j]*a[ixj];
        }
    }         /* ------- end of big loop ------- */

    ny  = n-1;
    it  = n*n;
    for (j=1; j<=ny; j++) {
        ia  = it - j;
        ib  = n  - j;
        ic  = n;
        for (k=1; k<=j; k++) {
            b[ib]  -= a[ia]*b[ic];
            ia     -= n;
            ic--;
        }
    }
}
/*--------------------------------------------------------------------*/
double fillin(double x, double *ab, double *ord, int no)
{
int i, m;
double y;

/*     *** FILLIN *** PARABOLIC INTERPOLATION
 *     FIND Y(X) FROM TABLE OF
 *     AB(N) AND ORD(N) CONTAINING NO POINTS.
 */
#define ANTRA(x1,x2,x3,x,y1,y2,y3)     (double)((y1)*((x) - (x2))*((x) - \
         (x3))/(((x1) - (x2))*((x1) - (x3))) + (y2)*((x) - (x1))*((x) - \
         (x3))/(((x2) - (x1))*((x2) - (x3))) + (y3)*((x) - (x1))*((x) - \
         (x2))/(((x3) - (x1))*((x3) - (x2))))

    if (x==ab[0]) {                     /* case: x on left edge */
        y = ord[0]; return y;
    }
    else {
        if (x<ab[0]) {                  /* case: x < left edge  */
            y = ANTRA(ab[0],ab[1],ab[2], x, ord[0],ord[1],ord[2]);
            return y;
        }
        else {              /**** case: x>ab[0] hereafter to the end  ****/
            if (x==ab[1]) {             /* case: x on ab[1]  */
                y = ord[1]; return y;
            }
            else {
                if (x<ab[1]) {          /* case: ab[0]<x<ab[1] */
                    y = ANTRA(ab[0],ab[1],ab[2], x, ord[0],ord[1],ord[2]);
                    return y;
                }
                else {      /**** case: x>ab[1]   ****/
                    for( i = 3; i <= no; i++ ){
                        m = i;
                        if (x==ab[i-1]) {     /* case: x on ab[i-1]  */
                            y = ord[i - 1]; return y;
                        }
                        else {
                            if (x<ab[i-1]) {  /* case: x<ab[i-1]  */
                                y = ANTRA(ab[m-3],ab[m-2],ab[m-1], x,
                                          ord[m -3],ord[m-2],ord[m-1]);
                                return y;
                            }
                        }
                    }
                    /* extrapolate beyond right edge */
                    y = ANTRA(ab[m-3],ab[m-2],ab[m-1], x,
                              ord[m -3],ord[m-2],ord[m-1]);
                    return y;
                }
            }
        }
    }
#undef  ANTRA
}
/*--------------------------------------------------------------------*/
/* powi.c - calc x^n, where n is an integer! */

#define IS_ODD(j)       ((j) & 1 )

double powi( double x, int n )         /* returns:  x^n */
        /* x: base,  n: exponent */
{
    double p;               /* holds partial product */

    if (x == 0)
        return( 0. );
    if (n < 0){           /* test for negative exponent */
        n = -n;
        x = 1/x;
    }

    p = IS_ODD(n) ? x : 1;  /* test & set zero power */

    while ( n >>= 1 ){       /* now do the other powers */
        x *= x;         /* sq previous power of x */
        if (IS_ODD(n)) /* if low order bit set */
                p *= x; /*then, multiply partial product by latest power of x*/
    }
    return (p);
}
#undef IS_ODD

//-----------------------------------------------------------------------
double powr( double x, double apower )         // returns:  x^apower 
        // x: base,  a: exponent 
{
    double p;               // holds partial product 

    if (x < 0.0) {
      cout << "powr: the base of powr is negative.\n";
      cout << "powr:   the base of powr is negative.\n";
      cout << "powr:     the base of powr is negative.\n";
      return 0.0;
    }
    if (x == 0.0) return( 0. );
    if (fabs(x-1.0)<0.00001f) {
      return 1.0;
    }
    p = apower * log(x);
    p = exp(p);
    return (p);
}
/*----------------------------------------------------------------------*/
/*-----Brute FORCE Version of 6/16/81 by D. Greeley---------------------*
 *     C-Version by LCS, March 1992.
 *----------------------------------------------------------------------*/
void vorseg(double xp, double yp, double zp, double  x1, double y1, double z1,
            double x2, double y2, double z2, double *vx, double *vy, double *vz,
            double *sx, double *sy, double *sz, int k)
{
    /*----Single precision for suburban/far_field----*/
    /*----Change to double if vorseg_length is longer than approx. 20---*/
    double   a, asq, ax, ay, az, b, bsq, c, csq, d, dsq, dx, dy, dz,
            e, eps, f, foasq, fsq, fx, fy, fz, hx, hy, hz, t, w;//..tol
    double   tol1, tol2, tol_preset=0.0004; /*--0.0025 if float--*/ /*--0.0004 if double--*/
    double  ffrsq=25.0, nfrsq=0.01, half=0.5, two=2.0, four=4.0;
    /*------Double precision variables for near_field case-----*/
    double ad, asqd, axd, ayd, azd, bd, bsqd, cd, csqd, dd, dsqd,
           dxd, dyd, dzd, ed, epsd, fd, fsqd, fxd, fyd, fzd, hxd,
           hyd, hzd, td, wd;
    double halfd=0.5, twod=2.0;

    /*-----VORTEX(A),MEAN DISTANCE(F),AND NORMAL(H) VECTORS------*/
    /*...nfrsq=(nearfield radius)^2, ffrsq=(farfield radius)^2...*/
    ax = x2 - x1;
    ay = y2 - y1;
    az = z2 - z1;
    asq = ax*ax + ay*ay + az*az;
    fx = half * (x2 + x1) - xp;
    fy = half * (y2 + y1) - yp;
    fz = half * (z2 + z1) - zp;
    fsq = fx*fx + fy*fy + fz*fz;
    foasq = fsq/asq;
    if (foasq < nfrsq && k == 1)
        goto near_field;  /*--Case when f.p. is near vorseg center--*/

    hx = az*fy - ay*fz;
    hy = ax*fz - az*fx;
    hz = ay*fx - ax*fy;
    if (foasq < ffrsq)
        goto suburban_field;

    /*-----FAR FIELD APPROXIMATIONS------ */
    f = sqrt(fsq);
    t = half/(f*f*f);
    *vx = t * hx;
    *vy = t * hy;
    *vz = t * hz;
    if (k==0)
            return;
    a = sqrt(asq);
    t *= -a;
    *sx = t * fx;
    *sy = t * fy;
    *sz = t * fz;
    return;

    /*-----NEAR AND SURBURBAN FIELD FOR VORTEX ONLY AND ---------------------
     *-----SUBURBAN FIELD FOR VORTEX AND SOURCE ----------------------------- */
suburban_field:
    bsq = (x2-xp)*(x2-xp) + (y2-yp)*(y2-yp) + (z2-zp)*(z2-zp);
    csq = (x1-xp)*(x1-xp) + (y1-yp)*(y1-yp) + (z1-zp)*(z1-zp);
    dsq = (hx*hx + hy*hy + hz*hz)/asq;
    a = sqrt(asq);
    b = sqrt(bsq);
    c = sqrt(csq);
    d = sqrt(dsq);
    e = (asq + csq - bsq)/(two*a);
    /*------check if f.p. falls on/extension of vorseg end points---*/
    if (dsq <= 0.0) {
        t = fabs( a*(a - 2.0*e)/(4.0*e*e*(a-e)*(a-e)) ) /a;
    }
    else {
        if (e >= 0.0 && e <= a) {
            t = ((a - e)/b + e/c)/(2.0*dsq) /a;
        }
        else {
            tol1 = fabs(e)*tol_preset;
            tol2 = fabs(e-a)*tol_preset;
            if ( ((e<0.0) && (dsq<tol1)) || ((e>0.0) && (dsq<tol2)) ) {
                t = fabs( a*(a - 2.0*e)/(4.0*e*e*(a-e)*(a-e)) ) /a;
            }
            else {
                t = ((a - e)/b + e/c)/(2.0*dsq) /a;
            }
        }
    }

    *vx = t*hx;
    *vy = t*hy;
    *vz = t*hz;
    if (k == 0)
        return;
    t *= a;
    w = (c - b)/(two*a*b*c);
    eps = e/a;
    dx = x1 + eps*ax - xp;
    dy = y1 + eps*ay - yp;
    dz = z1 + eps*az - zp;
    *sx = w*ax - t*dx;
    *sy = w*ay - t*dy;
    *sz = w*az - t*dz;
    return;

    /*-----NEAR FIELD (VORTEX & SOURCE) ------*/
    /*----- Double Precision Calculation -----*/
near_field:
    axd = x2 - x1;
    ayd = y2 - y1;
    azd = z2 - z1;
    asqd = axd*axd + ayd*ayd + azd*azd;
    ad = sqrt(asqd);
    fxd = halfd*(x2 + x1) - xp;
    fyd = halfd*(y2 + y1) - yp;
    fzd = halfd*(z2 + z1) - zp;
    fsqd = fxd*fxd + fyd*fyd + fzd*fzd;
    fd = sqrt(fsqd);
    hxd = azd*fyd - ayd*fzd;
    hyd = axd*fzd - azd*fxd;
    hzd = ayd*fxd - axd*fyd;
    bsqd = (x2-xp)*(x2-xp) + (y2-yp)*(y2-yp) + (z2-zp)*(z2-zp);
    csqd = (x1-xp)*(x1-xp) + (y1-yp)*(y1-yp) + (z1-zp)*(z1-zp);
    dsqd = (hxd*hxd + hyd*hyd + hzd*hzd)/asqd;
    bd = sqrt(bsqd);
    cd = sqrt(csqd);
    dd = sqrt(dsqd);
    ed = (asqd + csqd - bsqd)/(twod*ad);
    td = ((ad - ed)/bd + ed/cd)/(twod*ad*dsqd);
    *vx = td * hxd;
    *vy = td * hyd;
    *vz = td * hzd;
    td *= ad;
    wd = (cd - bd)/(twod*ad*bd*cd);
    epsd = ed/ad;
    dxd = x1 + epsd*axd - xp;
    dyd = y1 + epsd*ayd - yp;
    dzd = z1 + epsd*azd - zp;
    *sx = wd*axd - td*dxd;
    *sy = wd*ayd - td*dyd;
    *sz = wd*azd - td*dzd;
    return;
}
/*----------------------------------------------------------------------*/
void prcros(double a[], double b[], double c[])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return;
}

double prdot(double a[], double b[])
{
    double prdot_v;

    prdot_v = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return( prdot_v );
}

void prcrosK(double *a, double *b, double *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return;
}

double prdotK(double *a, double *b)
{
    double prdot_v;

    prdot_v = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return( prdot_v );
}

/*----------------------------------------------------------------------*/
double potder(double xp, double xin[], double yin[])
{
        double potder_v, x1, x2, a1, b1;

        x1 = (yin[1] - yin[0])/(xin[1] - xin[0]);
        x2 = (yin[1] - yin[2])/(xin[1] - xin[2]);
        a1 = (x1 - x2)/(xin[0] - xin[2]);
        b1 = (yin[1] - yin[2])/(xin[1] - xin[2]) - a1*(xin[1] + xin[2]);
        potder_v = 2.0*a1* xp + b1;
        return( potder_v );
}

/*----------------------------------------------------------------------*/

/* ----------------------------------------------------------------------
 *     EXPSPL.c : COLLECTION OF THE SUBROUTINES ASSOCIATED WITH THE
 *                EXPONENTIAL(TENSION) SPLINES
 *
 *
 *     SUBROUTINE SPL PROVIDES THE CURVE FITTING BY USING EXPONENTIAL
 *     CUBIC SPLINE (SO-CALLED TENSION SPLINE) BASED ON McCARTIN'S PAPER.
 *     (McCARTIN B. J. "APPLICATION OF EXPONENTIAL SPLINE IN COMPUTATIONAL
 *      FLUID DYNAMICS", AIAA JOURNAL, 1983)
 *
 *     PROGRAMMED BY JUNG-CHUN SUH AND MYUNG-HWAN KIM, 1987
 *
 *...........................................................................
 *     INPUT ARGUMENT DESCRIPTION
 *       npt:   NO. OF POINTS TO BE FITTED (i.e. NO. OF GIVEN NODAL POINTS)
 *              NOT GREATER THAN NDIM (=146)
 *       xin[]: ABSCISSAS OF NODAL POINTS (NPT ASCENDING ARRAYED VALUES)
 *       yin[]: ORDINATES OF NODAL POINTS (NPT CORRESPONDING ARRAYED VALUES)
 *       slpl:  SLOPE AT 1ST (LEFT) POINT
 *              (BUILT-IN BY QUADRATIC FITTING OF ADJACENT POINTS)
 *       slpr:  SLOPE AT LAST (RIGHT) POINT
 *              (BUILT-IN BY QUADRATIC FITTING OF ADJACENT POINTS)
 *       relax: RELAXATION PARAMETER OF RELATIVE PORTION OF NEW TENSION VALUE
 *              TO PREVIOUS VALUE WHEN UPDATING TENSION PARAMETER
 *              (BETWEEN 0 AND 1) (RECOMMENDED VALUE: 0.5)
 *       itmax: ALLOWED NO. OF ITERATION DURING FINDING TENSION
 *              (RECOMMENDED VALUE: 20.)
 *       ispl:  CONTROL NO. FOR CHOICE OF SPLINE
 *              (=1: EXPONENTIAL CUBIC SPLINE, =0: CUBIC SPLINE)
 *...........................................................................
 *     OUTPUT ARGUMENTS
 *       r[]:   2ND DERIVATIVES AT NODAL POINTS
 *       p[]:   TENSION PARAMETERS APPLIED TO SOME INTEVALS
 *       *iter: NO. OF ITERATION TO BE CARRIED OUT TO FIND TENSION PARAMETER
 *              WITHIN EITHER GIVEN ROUND-OFF ALLOWANCE OR NO. OF THE ALLOWED
 *              ITERATIONS.
 *
 *     REQUIRED SUBROUTINES: TRIDAG (FOR SOLVING TRIDIAGONAL MATRIX SYSTEM)
 *     APPLICATION SUBROUTINES: ESPL (EVALUATION OF VALUES AT SOME POINTS)
 *                              DSPL ( DERIVATIVE OF FITTING CURVE)
 *                              INTSPL (INTEGRATION OF FITTING CURVE)
 *...........................................................................
 *     C-Version 1.1, March 18, 1993, LCS
 *...........................................................................
 */
#define NDIM    1000

void spl(int npt, double Xin[], double Yin[], double *slpl, double *slpr,
         double relax, int itmax, int ispl, double R[], double P[], int *iter)
{
    int i, i1, ich[NDIM+1], idb, it, iupdt, k, k1, k2,
        kep[NDIM+1], kk, n, nch, neip, nin;
    double a[NDIM+1], a1, a2, aa, alamda, anum, b[NDIM+1], bb, c[NDIM+1],
     d[NDIM+1], deno, e[NDIM+1], h[NDIM+1], ph, phmax, pn, po[NDIM+1], pp,
     rhs[NDIM+1], roundoff, sh;
    /* OFFSET Vectors w/subscript range: 1 to dimension */
    double *const xin = &Xin[0] - 1;
    double *const yin = &Yin[0] - 1;
    double *const r = &R[0] - 1;
    double *const p = &P[0] - 1;
    double *const A = &a[0] + 1;
    double *const B = &b[0] + 1;
    double *const C = &c[0] + 1;

    /*----Implicit functions for coefficients of
          quadratic fitting of 3 points. */
#define FF1(x1,x2,x3,y1,y2,y3) \
        ((((y1)-(y2))/((x1)-(x2))-((y2)-(y3))/((x2)-(x3)))/((x1)-(x3)))
#define FF2(x1,x2,y1,y2,aa)   \
           (((y1)-(y2))/((x1)-(x2))-(aa)*((x1)+(x2)))

    /*-----Slopes at end points-----*/
    aa = FF1( xin[1], xin[2], xin[3], yin[1], yin[2], yin[3] );
    bb = FF2( xin[1], xin[2], yin[1], yin[2], aa );
    *slpl = 2*aa*xin[1] + bb;
    aa = FF1( xin[npt-2], xin[npt-1], xin[npt], yin[npt-2],
              yin[npt-1], yin[npt] );
    bb = FF2( xin[npt-2], xin[npt-1], yin[npt-2], yin[npt-1],aa );
    *slpr = 2*aa*xin[npt] + bb;
    /*-----------------------------*/

    phmax = 50.;
    roundoff = 1.e-8;
    /* delta x & delta y */
    nin = npt - 1;
    for (i = 1; i <= nin; i++) {
        h[i] = xin[i + 1] - xin[i];
    }
    /* Tridiagonal matrix */
    b[1] = h[1]/3.;
    c[1] = h[1]/6.;
    r[1] = (yin[2] - yin[1])/h[1] - *slpl;
    for (n = 2; n <= nin; n++) {
        a[n] = h[n - 1]/6.;
        b[n] = (h[n - 1] + h[n])/3.;
        c[n] = h[n]/6.;
        r[n] = (yin[n+1] - yin[n])/h[n] - (yin[n] - yin[n-1])/h[n-1];
    }
    a[npt] = c[nin];
    b[npt] = h[nin]/3.;
    r[npt] = *slpr - (yin[npt] - yin[nin])/h[nin];
    /* Save RHS */
    for (n = 1; n <= npt; n++) {
        rhs[n] = r[n];
    }
    /*  Use cubic spline w/o tension as zeroth iteration */

    tridag(1, npt, A, B, C, R);
    *iter = 0;
    if (ispl == 0)
        return;

    /*     Exponential spline
     *     Check intervals where EIPs might be occured */
    nch = 0;
    for (i = 1; i <= nin; i++) {
        if (rhs[i]*rhs[i + 1] >= 0.) {
            nch += 1;
            ich[nch] = i;
        }
    }
    if (nch == 0)
        return;
    /* Find tension parameters by iteration */
    for (i = 1; i <= nin; i++) {
        d[i] = h[i]/3.;
        e[i] = h[i]/6.;
        p[i] = 0.;
    }

    for (it = 1; it <= itmax; it++) {
        k2 = 0;
        /*  printf(" \n");
            for (nn=1; nn<=npt; nn++) {
            printf(" %12.5e %12.5e %12.5e %12.5e \n",
                     xin[NN],rhs[NN],r[NN],p[NN]);
            }
        */
        iupdt = 0;
        /* Check extraneous inflection points */
        neip = 0;
        for (n = 1; n <= nch; n++) {
            i1 = ich[n];
            for (i = i1; i <= (i1 + 1); i++) {
                if ((r[i]*rhs[i] > 0. || (fabs( r[i] ) <= roundoff &&
                 fabs( rhs[i] ) <= roundoff)) || (neip >= 1 && i ==
                 kep[neip]))
                    goto L_60;
                neip += 1;
                kep[neip] = i;
L_60:
                ;
            }
        }
        if (neip == 0 && it == 1) {
            /* printf("..Coefficients are same as ones in cubic splines\n);*/
            *iter = it - 1;
            return;
        }
        if (neip == 0) {
            *iter = it - 1;
            return;
        }

        for (i = 1; i <= nin; i++) {
            po[i] = p[i];
        }
        for (n = 1; n <= neip; n++) {
            k = kep[n];
            /* Set lambda */
            a1 = fabs( rhs[k] );
            a2 = b[k]*fabs( r[k] );
            anum = fmax( a1, a2 );
            if (k == 1)
                deno = fabs( r[k + 1] );
            if (k == npt)
                deno = fabs( r[k - 1] );
            if (k > 1 && k < npt)
                deno = 2.*fmax( fabs( r[k - 1] ), fabs( r[k + 1] ) );
            if (deno <= roundoff || anum <= roundoff)
                goto L_90;
            alamda = anum/deno;
            /*  Update tension parameter */
            idb = 0;
            k1 = k - 1;
            if (k == 1)
                k1 = 1;
            if (n > 1 && k1 == k2)
                idb = 1;
            k2 = k;
            if (k == npt)
                k2 = nin;
            for (kk = k1; kk <= k2; kk++) {
                a1 = 1./sqrt( alamda*h[kk] );
                if (a1 <= po[kk])
                    goto L_90;
                /*  DP=relax*(a1-po[kk]) */
                pn = po[kk] + relax*(a1 - po[kk]);
                if (idb == 1 && kk == k1) {
                    if (pn < p[kk])
                        p[kk] = pn;
                } else{
                    p[kk] = pn;
                }
                /* Compute coefficients of system of eq. for 2nd derivatives */
                ph = p[kk]*h[kk];
                pp = p[kk]*p[kk];
                if (ph < phmax) {
                    sh = sinh( ph );
                    d[kk] = (p[kk]*cosh( ph )/sh - 1./h[kk])/pp;
                    e[kk] = (1./h[kk] - p[kk]/sh)/pp;
                } else{
                    d[kk] = (p[kk] - 1./h[kk])/pp;
                    e[kk] = (1./h[kk])/pp;
                }
            }
            iupdt = 1;
L_90:
            ;
        }
        if (iupdt == 0) {
            *iter = it - 1;
            return;
        }
        /* Update elements of tridiagonal matrix */
        b[1] = d[1];
        c[1] = e[1];
        for (i = 2; i <= nin; i++) {
            a[i] = e[i - 1];
            b[i] = d[i - 1] + d[i];
            c[i] = e[i];
        }
        a[npt] = e[nin];
        b[npt] = d[nin];
        for (n = 1; n <= npt; n++) {
            r[n] = rhs[n];
        }

        tridag(1, npt, A, B, C, R);

    }

    /*  printf(" No convergence in finding tension parameters\n");
        printf(" during %3d iterations\n", itmax);      */
    *iter = itmax;

    return;
#undef  FF2
#undef  FF1
}

void splK(int ispl, double relax, double slpl, double slpr, 
          int npt, double *Xin, double *Yin, double *R, double *P)
{
    int i, i1, ich[NDIM+1], idb, it, iupdt, k, k1, k2,
        kep[NDIM+1], kk, n, nch, neip, nin;
    double a[NDIM+1], a1, a2, aa, alamda, anum, b[NDIM+1], bb, c[NDIM+1],
     d[NDIM+1], deno, e[NDIM+1], h[NDIM+1], ph, phmax, pn, po[NDIM+1], pp,
     rhs[NDIM+1], roundoff, sh;
    ///////////////////////////
    int itmax = 20, iter;
    ///////////////////////////
    /* OFFSET Vectors w/subscript range: 1 to dimension */
    double *xin, *yin, *r, *p;
    xin = Xin - 1;
    yin = Yin - 1;
    r = R - 1;
    p = P - 1;
    double *const A = &a[0] + 1;
    double *const B = &b[0] + 1;
    double *const C = &c[0] + 1;

    /*----Implicit functions for coefficients of
          quadratic fitting of 3 points. */
#define FF1(x1,x2,x3,y1,y2,y3) \
        ((((y1)-(y2))/((x1)-(x2))-((y2)-(y3))/((x2)-(x3)))/((x1)-(x3)))
#define FF2(x1,x2,y1,y2,aa)   \
           (((y1)-(y2))/((x1)-(x2))-(aa)*((x1)+(x2)))

    /*-----Slopes at end points-----*/
    if (slpl == 0.0) {
       aa = FF1( xin[1], xin[2], xin[3], yin[1], yin[2], yin[3] );
       bb = FF2( xin[1], xin[2], yin[1], yin[2], aa );
       slpl = 2*aa*xin[1] + bb;
    }
    if (slpr == 0.0) {
       aa = FF1( xin[npt-2], xin[npt-1], xin[npt], yin[npt-2], yin[npt-1], yin[npt] );
       bb = FF2( xin[npt-2], xin[npt-1], yin[npt-2], yin[npt-1],aa );
       slpr = 2*aa*xin[npt] + bb;
    }
    /*-----------------------------*/

    phmax = 50.;
    roundoff = 1.e-8;
    /* delta x & delta y */
    nin = npt - 1;
    for (i = 1; i <= nin; i++) {
        h[i] = xin[i + 1] - xin[i];
    }
    /* Tridiagonal matrix */
    b[1] = h[1]/3.;
    c[1] = h[1]/6.;
    r[1] = (yin[2] - yin[1])/h[1] - slpl;
    for (n = 2; n <= nin; n++) {
        a[n] = h[n - 1]/6.;
        b[n] = (h[n - 1] + h[n])/3.;
        c[n] = h[n]/6.;
        r[n] = (yin[n+1] - yin[n])/h[n] - (yin[n] - yin[n-1])/h[n-1];
    }
    a[npt] = c[nin];
    b[npt] = h[nin]/3.;
    r[npt] = slpr - (yin[npt] - yin[nin])/h[nin];
    /* Save RHS */
    for (n = 1; n <= npt; n++) {
        rhs[n] = r[n];
    }
    /*  Use cubic spline w/o tension as zeroth iteration */

    tridag(1, npt, A, B, C, R);
    iter = 0;
    if (ispl == 0)
        return;

    /*     Exponential spline
     *     Check intervals where EIPs might be occured */
    nch = 0;
    for (i = 1; i <= nin; i++) {
        if (rhs[i]*rhs[i + 1] >= 0.) {
            nch += 1;
            ich[nch] = i;
        }
    }
    if (nch == 0)
        return;
    /* Find tension parameters by iteration */
    for (i = 1; i <= nin; i++) {
        d[i] = h[i]/3.;
        e[i] = h[i]/6.;
        p[i] = 0.;
    }

    for (it = 1; it <= itmax; it++) {
        k2 = 0;
        /*  printf(" \n");
            for (nn=1; nn<=npt; nn++) {
            printf(" %12.5e %12.5e %12.5e %12.5e \n",
                     xin[NN],rhs[NN],r[NN],p[NN]);
            }
        */
        iupdt = 0;
        /* Check extraneous inflection points */
        neip = 0;
        for (n = 1; n <= nch; n++) {
            i1 = ich[n];
            for (i = i1; i <= (i1 + 1); i++) {
                if ((r[i]*rhs[i] > 0. || (fabs( r[i] ) <= roundoff &&
                 fabs( rhs[i] ) <= roundoff)) || (neip >= 1 && i ==
                 kep[neip]))
                    goto L_60;
                neip += 1;
                kep[neip] = i;
L_60:
                ;
            }
        }
        if (neip == 0 && it == 1) {
            /* printf("..Coefficients are same as ones in cubic splines\n);*/
            iter = it - 1;
            return;
        }
        if (neip == 0) {
            iter = it - 1;
            return;
        }

        for (i = 1; i <= nin; i++) {
            po[i] = p[i];
        }
        for (n = 1; n <= neip; n++) {
            k = kep[n];
            /* Set lambda */
            a1 = fabs( rhs[k] );
            a2 = b[k]*fabs( r[k] );
            anum = fmax( a1, a2 );
            if (k == 1)
                deno = fabs( r[k + 1] );
            if (k == npt)
                deno = fabs( r[k - 1] );
            if (k > 1 && k < npt)
                deno = 2.*fmax( fabs( r[k - 1] ), fabs( r[k + 1] ) );
            if (deno <= roundoff || anum <= roundoff)
                goto L_90;
            alamda = anum/deno;
            /*  Update tension parameter */
            idb = 0;
            k1 = k - 1;
            if (k == 1)
                k1 = 1;
            if (n > 1 && k1 == k2)
                idb = 1;
            k2 = k;
            if (k == npt)
                k2 = nin;
            for (kk = k1; kk <= k2; kk++) {
                a1 = 1./sqrt( alamda*h[kk] );
                if (a1 <= po[kk])
                    goto L_90;
                /*  DP=relax*(a1-po[kk]) */
                pn = po[kk] + relax*(a1 - po[kk]);
                if (idb == 1 && kk == k1) {
                    if (pn < p[kk])
                        p[kk] = pn;
                } else{
                    p[kk] = pn;
                }
                /* Compute coefficients of system of eq. for 2nd derivatives */
                ph = p[kk]*h[kk];
                pp = p[kk]*p[kk];
                if (ph < phmax) {
                    sh = sinh( ph );
                    d[kk] = (p[kk]*cosh( ph )/sh - 1./h[kk])/pp;
                    e[kk] = (1./h[kk] - p[kk]/sh)/pp;
                } else{
                    d[kk] = (p[kk] - 1./h[kk])/pp;
                    e[kk] = (1./h[kk])/pp;
                }
            }
            iupdt = 1;
L_90:
            ;
        }
        if (iupdt == 0) {
            iter = it - 1;
            return;
        }
        /* Update elements of tridiagonal matrix */
        b[1] = d[1];
        c[1] = e[1];
        for (i = 2; i <= nin; i++) {
            a[i] = e[i - 1];
            b[i] = d[i - 1] + d[i];
            c[i] = e[i];
        }
        a[npt] = e[nin];
        b[npt] = d[nin];
        for (n = 1; n <= npt; n++) {
            r[n] = rhs[n];
        }

        tridag(1, npt, A, B, C, R);

    }

    /*  printf(" No convergence in finding tension parameters\n");
        printf(" during %3d iterations\n", itmax);      */
    iter = itmax;

    return;
#undef  FF2
#undef  FF1
}

/* ----------------------------------------------------------------------
 *     SUBROUTINE ESPL PROVIDES THE EVALUATION OF VALUES AT SOME POINTS
 *     FROM CURVE FITTED BY USING EXPONENTIAL(TENSION) CUBIC SPLINE
 *     BASED ON McCARTIN'S PAPER. THUS SUBROUTINE SPL SHOULD BE CALLED
 *     FIRST IN ORDER TO USE THIS SUBROUTINE.
 *     (McCARTIN B. J. "APPLICATION OF EXPONENTIAL SPLINE IN
 *        COMPUTATIONAL FLUID DYNAMICS", AIAA JOURNAL, 1983)
 *
 *     PROGRAMMED BY JUNG-CHUN SUH AND MYUNG-HWAN KIM, 1987
 *.........................................................................
 *     INPUT ARGUMENT DESCRIPTION
 *       npt:   NO. OF POINTS TO BE FITTED (I.E. NO. OF GIVEN NODAL POINTS)
 *       xin[]: ABSCISSAS OF NODAL POINTS (NPT ASCENDING ARRAYED VALUES)
 *       yin[]: ORDINATES OF NODAL POINTS (NPT CORRESPONDING ARRAYED VALUES)
 *       ispl:  CONTROL NO. FOR CHOICE OF SPLINE
 *              (=1: EXPONENTIAL CUBIC SPLINE, =0: CUBIC SPLINE)
 *       f2[]:  2ND DERIVATIVES AT NODAL POINTS, COMPUTED IN SUBROUTINE SPL
 *       p[]:   TENSION PARAMETERS APPLIED TO SOME INTEVALS, COMPUTED IN
 *              SUBROUTINE SPL
 *.........................................................................
 *     OUTPUT ARGEMENTS
 *       nout:   NO. OF POINTS TO BE EVALUATED
 *       xout[]: ABSCISSAS OF POINTS TO BE EVALUATED
 *       yout[]: EVALUATED ORDINATES OF POINTS
 *
 *     RELATED ROUTINES: SPL (FOR FITTING GIVEN DATA) WHICH SHOULD BE
 *                          CALLED BEFORE ESPL IS CALLED.
 *
 * ----------------------------------------------------------------------
 */
void espl(int npt, double Xin[], double Yin[], int ispl, double F2[],
          double P[], int nout, double Xout[], double Yout[])
{
    int m, n, nin, np1;
    double dx1, dx2, fa, fb, fb1, fb2, h, ph, phmax, pp, px1, px2,
          roundoff, xo;
    /* OFFSET Vectors w/subscript range: 1 to dimension */
    double *const xin  = &Xin[0] - 1;
    double *const yin  = &Yin[0] - 1;
    double *const f2   = &F2[0] - 1;
    double *const p    = &P[0] - 1;
    double *const xout = &Xout[0] - 1;
    double *const yout = &Yout[0] - 1;

    /*...............*/
    phmax = 50.;
    roundoff = 1.e-8;
    nin = npt - 1;
    for (m = 1; m <= nout; m++) {
        xo = xout[m];
        for (n = 1; n <= nin; n++) {
            if (xo <= xin[n + 1])
                goto L_20;
        }
        n = nin;
L_20:
        np1 = n + 1;
        dx1 = xo - xin[n];
        dx2 = xin[np1] - xo;
        h = dx1 + dx2;
        if (ispl == 1 && p[n] > roundoff) {
            ph = p[n]*h;
            px1 = p[n]*dx1;
            px2 = p[n]*dx2;
            pp = p[n]*p[n];
            fa = ((yin[n]-f2[n]/pp)*dx2 + (yin[np1]-f2[np1]/pp)*dx1)/h;
            if ((ph>phmax || fabs(px1)>phmax) || fabs(px2)>phmax) {
                fb1 = 0.;
                fb2 = 0.;
                if (fabs( fabs(px2)-ph ) < phmax)
                    fb1 = f2[n]*exp( fabs( px2 ) - ph );
                if (fabs( fabs(px1)-ph ) < phmax)
                    fb2 = f2[np1]*exp( fabs(px1) - ph );
                if (dx1 < 0.0)
                    fb2 = -fb2;
                if (dx2 < 0.0)
                    fb1 = -fb1;
                fb = (fb1 + fb2)/pp;
            } else{
                fb = (f2[n]*sinh(px2) + f2[np1]*sinh(px1))/(pp*sinh(ph));
            }
            yout[m] = fa + fb;
        } else{
            yout[m] = yin[n]*dx2/h + yin[np1]*dx1/h - f2[n]*dx2*(h -
                        dx2*dx2/h)/6. - f2[np1]*dx1*(h - dx1*dx1/h)/6.;
        }
    }
    return;
}

void esplK(int ispl, int npt, double *Xin, double *Yin, 
           int nout, double *Xout, double *Yout, double *F2, double *P)
{
    int m, n, nin, np1;
    double dx1, dx2, fa, fb, fb1, fb2, h, ph, phmax, pp, px1, px2,
          roundoff, xo;
    /* OFFSET Vectors w/subscript range: 1 to dimension */
    double *xin, *yin, *f2, *p, *xout, *yout;
    xin  = Xin - 1;
    yin  = Yin - 1;
    f2   = F2 - 1;
    p    = P - 1;
    xout = Xout - 1;
    yout = Yout - 1;

    /*...............*/
    phmax = 50.;
    roundoff = 1.e-8;
    nin = npt - 1;
    for (m = 1; m <= nout; m++) {
        xo = xout[m];
        for (n = 1; n <= nin; n++) {
            if (xo <= xin[n + 1])
                goto L_20;
        }
        n = nin;
L_20:
        np1 = n + 1;
        dx1 = xo - xin[n];
        dx2 = xin[np1] - xo;
        h = dx1 + dx2;
        if (ispl == 1 && p[n] > roundoff) {
            ph = p[n]*h;
            px1 = p[n]*dx1;
            px2 = p[n]*dx2;
            pp = p[n]*p[n];
            fa = ((yin[n]-f2[n]/pp)*dx2 + (yin[np1]-f2[np1]/pp)*dx1)/h;
            if ((ph>phmax || fabs(px1)>phmax) || fabs(px2)>phmax) {
                fb1 = 0.;
                fb2 = 0.;
                if (fabs( fabs(px2)-ph ) < phmax)
                    fb1 = f2[n]*exp( fabs( px2 ) - ph );
                if (fabs( fabs(px1)-ph ) < phmax)
                    fb2 = f2[np1]*exp( fabs(px1) - ph );
                if (dx1 < 0.0)
                    fb2 = -fb2;
                if (dx2 < 0.0)
                    fb1 = -fb1;
                fb = (fb1 + fb2)/pp;
            } else{
                fb = (f2[n]*sinh(px2) + f2[np1]*sinh(px1))/(pp*sinh(ph));
            }
            yout[m] = fa + fb;
        } else{
            yout[m] = yin[n]*dx2/h + yin[np1]*dx1/h - f2[n]*dx2*(h -
                        dx2*dx2/h)/6. - f2[np1]*dx1*(h - dx1*dx1/h)/6.;
        }
    }
    return;
}

/*----------------------------------------------------------------------
 *     SUBROUTINE DSPL PROVIDES THE 1ST AND 2ND DERIVATIVES AT SOME
 *     POINTS FROM CURVE FITTED BY USING EXPONENTIAL CUBIC SPLINE
 *     (SO-CALLED TENSION SPLINE) BASED ON McCARTIN'S PAPER.
 *     THUS SUBROUTINE SPL SHOULD BE CALLED FIRST IN ORDER TO USE
 *     THIS SUBROUTINE.
 *     (McCARTIN B. J. "APPLICATION OF EXPONENTIAL SPLINE IN
 *       COMPUTATIONAL FLUID DYNAMICS", AIAA JOURNAL, 1983)
 *
 *     PROGRAMMED BY JUNG-CHUN SUH AND MYUNG-HWAN KIM, 1987
 *.......................................................................
 *     INPUT ARGUMENT DESCRIPTION
 *       NPT: NO. OF POINTS TO BE FITTED (I.E. NO. OF GIVEN NODAL POINTS)
 *       XIN: ABSCISSAS OF NODAL POINTS (NPT ASCENDING ARRAYED VALUES)
 *       YIN: ORDINATES OF NODAL POINTS (NPT CORRESPONDING ARRAYED VALUES)
 *       ISPL: CONTROL NO. FOR CHOICE OF SPLINE
 *             (=1: EXPONENTIAL CUBIC SPLINE, =0: CUBIC SPLINE)
 *       F2: 2ND DERIVATIVES AT NODAL POINTS, COMPUTED IN SUBROUTINE SPL
 *       P: TENSION PARAMETERS APPLIED TO SOME INTEVALS, COMPUTED IN
 *          SUBROUTINE SPL
 *
 *.........................................................................
 *    OUTPUT ARGUMENTS
 *       nout: NO. OF POINTS TO BE EVALUATED
 *       xout: ABSCISSAS OF POINTS TO BE EVALUATED
 *       y1out: 1ST DERIVATIVE OF FITTED CURVE AT CORRESPONDING POINTS
 *       y2out: 2ND DERIVATIVE OF FITTED CURVE AT CORRESPONDING POINTS
 *
 *     RELATED ROUTINES: SPL (FOR FITTING GIVEN DATA) WHICH SHOULD BE
 *                       CALLED BEFORE DSPL IS CALLED.
 *
 *------------------------------------------------------------------------
 */
void dspl(int npt, double Xin[], double Yin[], int ispl, double F2[],
          double P[],  int nout, double Xout[], double Y1out[], double Y2out[])
{
    int m, n, nin, np1;
    double dx1, dx2, fa, fb, fb1, fb2, h, ph, phmax, pp, px1, px2,
          roundoff, xo;
    /* OFFSET Vectors w/subscript range: 1 to dimension */
    double *const xin   = &Xin[0] - 1;
    double *const yin   = &Yin[0] - 1;
    double *const f2    = &F2[0] - 1;
    double *const p     = &P[0] - 1;
    double *const xout  = &Xout[0] - 1;
    double *const y1out = &Y1out[0] - 1;
    double *const y2out = &Y2out[0] - 1;

    /*..........................*/
    phmax = 50.;
    roundoff = 1.e-8;
    nin = npt - 1;
    for (m = 1; m <= nout; m++) {
        xo = xout[m];
        for (n = 1; n <= nin; n++) {
            if (xo <= xin[n + 1])
                goto L_20;
        }
        n = nin;
L_20:
        np1 = n + 1;
        dx1 = xo - xin[n];
        dx2 = xin[np1] - xo;
        h = dx1 + dx2;
        if (ispl == 1 && p[n] > roundoff) {
            ph = p[n]*h;
            px1 = p[n]*dx1;
            px2 = p[n]*dx2;
            pp = p[n]*p[n];
            fa = (yin[np1] - yin[n] + (f2[n] - f2[np1])/pp)/h;
            if ((ph>phmax || fabs(px1)>phmax) || fabs(px2)>phmax) {
                fb1 = 0.;
                fb2 = 0.;
                if (fabs( fabs( px2 ) - ph ) < phmax)
                    fb1 = f2[n]*exp( fabs( px2 ) - ph );
                if (fabs( fabs( px1 ) - ph ) < phmax)
                    fb2 = f2[np1]*exp( fabs( px1 ) - ph );
                fb = (-fb1 + fb2)/p[n];
                if (dx1 < 0.0)
                    fb2 = -fb2;
                if (dx2 < 0.0)
                    fb1 = -fb1;
                y2out[m] = fb1 + fb2;
            } else{
                fb = (-f2[n]*cosh(px2) + f2[np1]*cosh(px1))/(p[n]*sinh(ph));
                y2out[m] = (f2[n]*sinh(px2) + f2[np1]*sinh(px1))/sinh(ph);
            }
            y1out[m] = fa + fb;
        } else{
            y1out[m] = (yin[np1] - yin[n])/h + (f2[n] - f2[np1])*h/6.
                        + (f2[np1]*dx1*dx1 - f2[n]*dx2*dx2)/(2.*h);
            y2out[m] = (f2[n]*dx2 + f2[np1]*dx1)/h;
        }
    }
    return;
}

void dsplK(int ispl, int npt, double *Xin, double *Yin,
           int nout, double *Xout, double *Y1out, double *Y2out, double *F2, double *P)
{
    int m, n, nin, np1;
    double dx1, dx2, fa, fb, fb1, fb2, h, ph, phmax, pp, px1, px2,
          roundoff, xo;
    /* OFFSET Vectors w/subscript range: 1 to dimension */
    double *xin, *yin, *f2, *p, *xout, *y1out, *y2out;
    xin   = Xin - 1;
    yin   = Yin - 1;
    f2    = F2 - 1;
    p     = P - 1;
    xout  = Xout - 1;
    y1out = Y1out - 1;
    y2out = Y2out - 1;

    /*..........................*/
    phmax = 50.;
    roundoff = 1.e-8;
    nin = npt - 1;
    for (m = 1; m <= nout; m++) {
        xo = xout[m];
        for (n = 1; n <= nin; n++) {
            if (xo <= xin[n+1])
                goto L_20;
        }
        n = nin;
L_20:
        np1 = n + 1;
        dx1 = xo - xin[n];
        dx2 = xin[np1] - xo;
        h = dx1 + dx2;
        if (ispl == 1 && p[n] > roundoff) {
            ph = p[n]*h;
            px1 = p[n]*dx1;
            px2 = p[n]*dx2;
            pp = p[n]*p[n];
            fa = (yin[np1] - yin[n] + (f2[n] - f2[np1])/pp)/h;
            if ((ph>phmax || fabs(px1)>phmax) || fabs(px2)>phmax) {
                fb1 = 0.;
                fb2 = 0.;
                if (fabs( fabs( px2 ) - ph ) < phmax)
                    fb1 = f2[n]*exp( fabs( px2 ) - ph );
                if (fabs( fabs( px1 ) - ph ) < phmax)
                    fb2 = f2[np1]*exp( fabs( px1 ) - ph );
                fb = (-fb1 + fb2)/p[n];
                if (dx1 < 0.0)
                    fb2 = -fb2;
                if (dx2 < 0.0)
                    fb1 = -fb1;
                y2out[m] = fb1 + fb2;
            } else{
                fb = (-f2[n]*cosh(px2) + f2[np1]*cosh(px1))/(p[n]*sinh(ph));
                y2out[m] = (f2[n]*sinh(px2) + f2[np1]*sinh(px1))/sinh(ph);
            }
            y1out[m] = fa + fb;
        } else{
            y1out[m] = (yin[np1] - yin[n])/h + (f2[n] - f2[np1])*h/6.
                        + (f2[np1]*dx1*dx1 - f2[n]*dx2*dx2)/(2.*h);
            y2out[m] = (f2[n]*dx2 + f2[np1]*dx1)/h;
        }
    }
    return;
}

/*------------------------------------------------------------------------
 *     SUBROUTINE INTSPL PROVIDES THE INTEGRATION FOR GIVEN INTERVAL
 *     FROM CURVE FITTED BY USING EXPONENTIAL CUBIC SPLINE (SO-CALLED
 *     TENSION SPLINE) BASED ON McCARTIN'S PAPER.
 *     THUS SUBROUTINE SPL SHOULD BE CALLED FIRST IN ORDER TO USE
 *     THIS SUBROUTINE.
 *     (McCARTIN B. J. "APPLICATION OF EXPONENTIAL SPLINE IN
 *       COMPUTATIONAL FLUID DYNAMICS", AIAA JOURNAL, 1983)
 *
 *     PROGRAMMED BY JUNG-CHUN SUH AND MYUNG-HWAN KIM, 1988
 *.........................................................................
 *     INPUT ARGUMENT DESCRIPTION
 *       npt:   NO. OF POINTS TO BE FITTED (I.E. NO. OF GIVEN NODAL POINTS)
 *       xin[]: ABSCISSAS OF NODAL POINTS (NPT ASCENDING ARRAYED VALUES)
 *       yin[]: ORDINATES OF NODAL POINTS (NPT CORRESPONDING ARRAYED VALUES)
 *       ispl:  CONTROL NO. FOR CHOICE OF SPLINE
 *              (=1: EXPONENTIAL CUBIC SPLINE, =0: CUBIC SPLINE)
 *       f2[]:  2ND DERIVATIVES AT NODAL POINTS, COMPUTED IN SUBROUTINE SPL
 *       p[]:   TENSION PARAMETERS APPLIED TO SOME INTEVALS, COMPUTED IN
 *              function spl().
 *.........................................................................
 *    OUTPUT ARGUMENTS
 *       xa:   ONE LIMIT (LOWER LIMIT) OF INTEGRAL
 *       xb:   THE OTHER LIMIT (UPPER LIMIT) OF INTEGRAL
 *       yint: INTEGRAL VALUE
 *
 *     RELATED ROUTINES: SPL (FOR FITTING GIVEN DATA) WHICH SHOULD BE
 *                       CALLED BEFORE INTSPL IS CALLED.
 *------------------------------------------------------------------------
 */
void intspl(int npt, double Xin[], double Yin[], int ispl,
            double F2[], double P[], double xa, double xb, double *yint)
{
    int ia, iap1, ib, ibp1, isign, j, jab, k, kp1, nin;
    double aa, bb, dum1, dum2, dx1, dx2, fa1, fa2, fa3, fa4, fl3, fl4,
          fr3, fr4, h, ph, phmax, pp, ppp, px1, px2, roundoff,
          yaout, ybout, ylout, yrout;
    /* OFFSET Vectors w/subscript range: 1 to dimension */
    double *const f2 = &F2[0] - 1;
    double *const p = &P[0] - 1;
    double *const xin = &Xin[0] - 1;
    double *const yin = &Yin[0] - 1;

    /*.....................*/
    phmax = 50.;
    roundoff = 1.e-8;
    if (xa == xb) {
        *yint = 0.;
        return;
    } else{
        if (xb > xa) {
            aa = xa;
            bb = xb;
            isign = 1;
        } else{
            aa = xb;
            bb = xa;
            isign = -1;
        }
    }

    nin = npt - 1;

    for (ia = 1; ia <= nin; ia++) {
        if (aa <= xin[ia + 1])
            goto L_20;
    }
    ia = nin;
L_20:
    iap1 = ia + 1;
    for (ib = 1; ib <= nin; ib++) {
        if (bb <= xin[ib + 1])
            goto L_22;
    }
    ib = nin;
L_22:
    ibp1 = ib + 1;

    jab = ib - ia;

    dx1 = aa - xin[ia];
    dx2 = xin[iap1] - aa;
    h = dx1 + dx2;
    if (ispl == 1 && p[ia] > roundoff) {
        ph = p[ia]*h;
        px1 = p[ia]*dx1;
        px2 = p[ia]*dx2;
        pp = p[ia]*p[ia];
        ppp = pp*p[ia];
        fa3 = (-(yin[ia] - f2[ia]/pp)*(dx2*dx2)
               + (yin[iap1] - f2[iap1]/pp)*(dx1*dx1))/(2.*h);
        fr3 = (yin[iap1] - f2[iap1]/pp)*h/2.;
        if ((ph>phmax || fabs(px1)>phmax) || fabs(px2)>phmax) {
            fa1 = 0.;
            fa2 = 0.;
            dum2 = fabs( px2 ) - ph - log( ppp );
            dum1 = fabs( px1 ) - ph - log( ppp );
            if (fabs( dum2 ) < phmax)
                fa1 = exp( dum2 );
            if (fabs( dum1 ) < phmax)
                fa2 = exp( dum1 );
            fa4 = -fa1*f2[ia] + fa2*f2[iap1];
        } else{
            fa4 = (-f2[ia]*cosh(px2) + f2[iap1]*cosh(px1))/(ppp*sinh(ph));
        }
        if (ph > phmax) {
            fr4 = f2[iap1]/ppp;
        } else{
            fr4 = (-f2[ia] + f2[iap1]*cosh( ph ))/(ppp*sinh( ph ));
        }
        yaout = fa4 + fa3;
        yrout = fr3 + fr4;
    } else{
        yaout = (-yin[ia]*(dx2*dx2) + yin[iap1]*(dx1*dx1))/(2.*h)
                + f2[ia]/24.*(dx2*dx2)*(2.*h - (dx2*dx2)/h)
                + f2[iap1]/24.*(dx1*dx1)*(-2.*h + (dx1*dx1)/h);
        yrout = yin[iap1]*h/2. - 1./24.*f2[iap1]*(h*h*h);
    }

    dx1 = bb - xin[ib];
    dx2 = xin[ibp1] - bb;
    h = dx1 + dx2;
    if (ispl == 1 && p[ib] > roundoff) {
        ph = p[ib]*h;
        px1 = p[ib]*dx1;
        px2 = p[ib]*dx2;
        pp = p[ib]*p[ib];
        ppp = pp*p[ib];
        fa3 = (-(yin[ib] - f2[ib]/pp)*(dx2*dx2)
               + (yin[ibp1] - f2[ibp1]/pp)*(dx1*dx1))/(2.*h);
        fl3 = -(yin[ib] - f2[ib]/pp)*h/2.;
        if ((ph > phmax || fabs( px1 ) > phmax) || fabs( px2 ) > phmax) {
            fa1 = 0.;
            fa2 = 0.;
            dum2 = fabs( px2 ) - ph - log( ppp );
            dum1 = fabs( px1 ) - ph - log( ppp );
            if (fabs( dum2 ) < phmax)
                fa1 = exp( dum2 );
            if (fabs( dum1 ) < phmax)
                fa2 = exp( dum1 );
            fa4 = -fa1*f2[ib] + fa2*f2[ibp1];
        } else{
            fa4 = (-f2[ib]*cosh(px2) + f2[ibp1]*cosh(px1))/(ppp*sinh(ph));
        }
        if (ph > phmax) {
            fl4 = -f2[ib]/ppp;
        } else{
            fl4 = (-f2[ib]*cosh(ph) + f2[ibp1])/(ppp*sinh(ph));
        }
        ybout = fa4 + fa3;
        ylout = fl3 + fl4;
    } else{
        ybout = (-yin[ib]*(dx2*dx2) + yin[ibp1]*(dx1*dx1))/(2.*h)
                 + f2[ib]/24.*(dx2*dx2)*(2.*h - (dx2*dx2)/h)
                 + f2[ibp1]/24.*(dx1*dx1)*(-2.*h + (dx1*dx1)/h);
        ylout = -yin[ib]*h/2. + 1./24.*f2[ib]*(h*h*h);
    }

    if (jab == 0)
        *yint = ybout - yaout;
    if (jab >= 1)
        *yint = yrout - yaout + ybout - ylout;
    if (jab > 1) {

        for (j = 1; j <= (jab - 1); j++) {
            k = ia + j;
            kp1 = k + 1;
            h = xin[kp1] - xin[k];
            if (ispl == 1 && p[k] > roundoff) {
                ph = p[k]*h;
                pp = p[k]*p[k];
                ppp = pp*p[k];
                fr3 = (yin[kp1] - f2[kp1]/pp)*h/2.;
                fl3 = -(yin[k] - f2[k]/pp)*h/2.;
                if (ph > phmax) {
                    fr4 = f2[kp1]/ppp;
                    fl4 = -f2[k]/ppp;
                } else{
                    fr4 = (-f2[k] + f2[kp1]*cosh( ph ))/(ppp*sinh( ph ));
                    fl4 = (-f2[k]*cosh( ph ) + f2[kp1])/(ppp*sinh( ph ));
                }
                yrout = fr3 + fr4;
                ylout = fl3 + fl4;
            } else{
                yrout = yin[kp1]*h/2. - 1./24.*f2[kp1]*(h*h*h);
                ylout = -yin[k]*h/2. + 1./24.*f2[k]*(h*h*h);
            }
            *yint += yrout - ylout;
        }

    }
    *yint *= isign;
    return;
}

void intsplK(int ispl, int npt, double *Xin, double *Yin, 
            double xa, double xb, double &yint, double *F2, double *P)
{
    int ia, iap1, ib, ibp1, isign, j, jab, k, kp1, nin;
    double aa, bb, dum1, dum2, dx1, dx2, fa1, fa2, fa3, fa4, fl3, fl4,
          fr3, fr4, h, ph, phmax, pp, ppp, px1, px2, roundoff,
          yaout, ybout, ylout, yrout;
    /* OFFSET Vectors w/subscript range: 1 to dimension */
    double *xin, *yin, *f2, *p;
    xin = Xin - 1;
    yin = Yin - 1;
    f2  = F2 - 1;
    p   = P - 1;


    /*.....................*/
    phmax = 50.;
    roundoff = 1.e-8;
    if (xa == xb) {
        yint = 0.;
        return;
    } else{
        if (xb > xa) {
            aa = xa;
            bb = xb;
            isign = 1;
        } else{
            aa = xb;
            bb = xa;
            isign = -1;
        }
    }

    nin = npt - 1;

    for (ia = 1; ia <= nin; ia++) {
        if (aa <= xin[ia + 1])
            goto L_20;
    }
    ia = nin;
L_20:
    iap1 = ia + 1;
    for (ib = 1; ib <= nin; ib++) {
        if (bb <= xin[ib + 1])
            goto L_22;
    }
    ib = nin;
L_22:
    ibp1 = ib + 1;

    jab = ib - ia;

    dx1 = aa - xin[ia];
    dx2 = xin[iap1] - aa;
    h = dx1 + dx2;
    if (ispl == 1 && p[ia] > roundoff) {
        ph = p[ia]*h;
        px1 = p[ia]*dx1;
        px2 = p[ia]*dx2;
        pp = p[ia]*p[ia];
        ppp = pp*p[ia];
        fa3 = (-(yin[ia] - f2[ia]/pp)*(dx2*dx2)
               + (yin[iap1] - f2[iap1]/pp)*(dx1*dx1))/(2.*h);
        fr3 = (yin[iap1] - f2[iap1]/pp)*h/2.;
        if ((ph>phmax || fabs(px1)>phmax) || fabs(px2)>phmax) {
            fa1 = 0.;
            fa2 = 0.;
            dum2 = fabs( px2 ) - ph - log( ppp );
            dum1 = fabs( px1 ) - ph - log( ppp );
            if (fabs( dum2 ) < phmax)
                fa1 = exp( dum2 );
            if (fabs( dum1 ) < phmax)
                fa2 = exp( dum1 );
            fa4 = -fa1*f2[ia] + fa2*f2[iap1];
        } else{
            fa4 = (-f2[ia]*cosh(px2) + f2[iap1]*cosh(px1))/(ppp*sinh(ph));
        }
        if (ph > phmax) {
            fr4 = f2[iap1]/ppp;
        } else{
            fr4 = (-f2[ia] + f2[iap1]*cosh( ph ))/(ppp*sinh( ph ));
        }
        yaout = fa4 + fa3;
        yrout = fr3 + fr4;
    } else{
        yaout = (-yin[ia]*(dx2*dx2) + yin[iap1]*(dx1*dx1))/(2.*h)
                + f2[ia]/24.*(dx2*dx2)*(2.*h - (dx2*dx2)/h)
                + f2[iap1]/24.*(dx1*dx1)*(-2.*h + (dx1*dx1)/h);
        yrout = yin[iap1]*h/2. - 1./24.*f2[iap1]*(h*h*h);
    }

    dx1 = bb - xin[ib];
    dx2 = xin[ibp1] - bb;
    h = dx1 + dx2;
    if (ispl == 1 && p[ib] > roundoff) {
        ph = p[ib]*h;
        px1 = p[ib]*dx1;
        px2 = p[ib]*dx2;
        pp = p[ib]*p[ib];
        ppp = pp*p[ib];
        fa3 = (-(yin[ib] - f2[ib]/pp)*(dx2*dx2)
               + (yin[ibp1] - f2[ibp1]/pp)*(dx1*dx1))/(2.*h);
        fl3 = -(yin[ib] - f2[ib]/pp)*h/2.;
        if ((ph > phmax || fabs( px1 ) > phmax) || fabs( px2 ) > phmax) {
            fa1 = 0.;
            fa2 = 0.;
            dum2 = fabs( px2 ) - ph - log( ppp );
            dum1 = fabs( px1 ) - ph - log( ppp );
            if (fabs( dum2 ) < phmax)
                fa1 = exp( dum2 );
            if (fabs( dum1 ) < phmax)
                fa2 = exp( dum1 );
            fa4 = -fa1*f2[ib] + fa2*f2[ibp1];
        } else{
            fa4 = (-f2[ib]*cosh(px2) + f2[ibp1]*cosh(px1))/(ppp*sinh(ph));
        }
        if (ph > phmax) {
            fl4 = -f2[ib]/ppp;
        } else{
            fl4 = (-f2[ib]*cosh(ph) + f2[ibp1])/(ppp*sinh(ph));
        }
        ybout = fa4 + fa3;
        ylout = fl3 + fl4;
    } else{
        ybout = (-yin[ib]*(dx2*dx2) + yin[ibp1]*(dx1*dx1))/(2.*h)
                 + f2[ib]/24.*(dx2*dx2)*(2.*h - (dx2*dx2)/h)
                 + f2[ibp1]/24.*(dx1*dx1)*(-2.*h + (dx1*dx1)/h);
        ylout = -yin[ib]*h/2. + 1./24.*f2[ib]*(h*h*h);
    }

    if (jab == 0)
        yint = ybout - yaout;
    if (jab >= 1)
        yint = yrout - yaout + ybout - ylout;
    if (jab > 1) {

        for (j = 1; j <= (jab - 1); j++) {
            k = ia + j;
            kp1 = k + 1;
            h = xin[kp1] - xin[k];
            if (ispl == 1 && p[k] > roundoff) {
                ph = p[k]*h;
                pp = p[k]*p[k];
                ppp = pp*p[k];
                fr3 = (yin[kp1] - f2[kp1]/pp)*h/2.;
                fl3 = -(yin[k] - f2[k]/pp)*h/2.;
                if (ph > phmax) {
                    fr4 = f2[kp1]/ppp;
                    fl4 = -f2[k]/ppp;
                } else{
                    fr4 = (-f2[k] + f2[kp1]*cosh( ph ))/(ppp*sinh( ph ));
                    fl4 = (-f2[k]*cosh( ph ) + f2[kp1])/(ppp*sinh( ph ));
                }
                yrout = fr3 + fr4;
                ylout = fl3 + fl4;
            } else{
                yrout = yin[kp1]*h/2. - 1./24.*f2[kp1]*(h*h*h);
                ylout = -yin[k]*h/2. + 1./24.*f2[k]*(h*h*h);
            }
            yint += yrout - ylout;
        }

    }
    yint *= isign;
    return;
}

/*-------------------------------------------------------------------*/
void tridag(int il, int iu, double A[], double B[], double C[], double R[])
{
    int i, j, lp;
    double bb[NDIM+1], q;
        /* OFFSET Vectors w/subscript range: 1 to dimension */
    double *const a = &A[0] - 1;
    double *const b = &B[0] - 1;
    double *const c = &C[0] - 1;
    double *const r = &R[0] - 1;

    /*  Tridiagonal matrix systme solver  */
    lp = il + 1;
    bb[il] = b[il];
    for (i = lp; i <= iu; i++) {
        q = a[i]/bb[i - 1];
        bb[i] = b[i] - q*c[i - 1];
        r[i] += -q*r[i - 1];
    }
    r[iu] /= bb[iu];
    for (i = lp; i <= iu; i++) {
        j = iu - i + il;
        r[j] = (r[j] - c[j]*r[j + 1])/bb[j];
    }
    return;
}
#undef NDIM                          

                                                      
//...........matrix memory allocation....12/28/96...LCS
int *ivec(int n_row) {
  int i;
  int *m;

  m = new int [n_row];
  if (!m) cout << "Insufficient memory in vector()\n";
  for (i=0; i<n_row; i++) {
    m[i] = 0;
  }
  return m;
}

double *vec(int n_row) {
  int i;
  double *m;

  m = new double [n_row];
  if (!m) cout << "Insufficient memory in vector()\n";
  for (i=0; i<n_row; i++) {
    m[i] = 0.0;
  }
  return m;
}

double **matrix2d(int n_row, int n_col) {
  int i, j;
  double **m;

  m = new double *[n_row];
  if (!m) cout << "Insufficient memory in matrix2d()\n";
  for (i=0; i<n_row; i++) {
    m[i] = new double [n_col];
    if (!m[i]) cout << "Insufficient memory in matrix2d()\n";
    for (j=0; j<n_col; j++) {
      m[i][j] = 0.0;
    }
  }
  return m;
}
void free_matrix2d(double **m, int n_row) {
  int i;
  for (i=0; i<n_row; i++) {
    delete m[i];
  }
  delete m;
}

//...integer matrix2d allocation...
int **imatrix2d(int n_row, int n_col) {
  int i, j;
  int **m;

  m = new int *[n_row];
  if (!m) cout << "Insufficient memory in int_matrix2d()\n";
  for (i=0; i<n_row; i++) {
    m[i] = new int[n_col];
    if (!m[i]) cout << "Insufficient memory in int_matrix2d()\n";
    for (j=0; j<n_col; j++) {
      m[i][j] = 0;
    }
  }
  return m;
}
void free_imatrix2d(int **m, int n_row) {
  int i;
  for (i=0; i<n_row; i++) {
    delete m[i];
  }
  delete m;
}

//...3D array allocation...
double ***matrix3d(int n_page, int n_row, int n_col) {
  int i, j, k;
  double ***m;
  m = new double **[n_page];
  if (!m) cout << "Insufficient memory in 3D_matrix3d()\n";
  for (i=0; i<n_page; i++) {
    m[i] = new double *[n_row];
    if (!m[i]) cout << "Insufficient memory in matrix3d()\n";
	  for (j=0; j<n_row; j++) {
      m[i][j] = new double [n_col];
      if (!m[i][j]) cout << "Insufficient memory in matrix3d()\n";
      for (k=0; k<n_col; k++) {
        m[i][j][k] = 0.0;
      }
    }
  }
  return m;
}

void free_matrix3d(double ***m, int n_page, int n_row) {
  int i, j;
  for (i=0; i<n_page; i++) {
    for (j=0; j<n_row; j++) {
      delete m[i][j];
    }
    delete m[i];
  }
  delete m;
}

//...4D array allocation...
double ****matrix4d(int n_box, int n_page, int n_row, int n_col) {
  int i, j, k, n;
  double ****m;
  m = new double ***[n_box];
  if (!m) cout << "Insufficient memory in matrix4d()\n";
  for (i=0; i<n_box; i++) {
    m[i] = new double **[n_page];
    if (!m[i]) cout << "Insufficient memory in matrix4d()\n";
	  for (j=0; j<n_page; j++) {
      m[i][j] = new double *[n_row];
      if (!m[i][j]) cout << "Insufficient memory in matrix3d()\n";
      for (k=0; k<n_row; k++) {
        m[i][j][k] = new double[n_col];
				if (!m[i][j][k]) cout << "Insufficient memory in matrix4d()\n";
        for (n=0; n<n_col; n++) {
          m[i][j][k][n] = 0.0;
        }
      }//...k
    }//...j
  }//...i
  return m;
}

void free_matrix4d(double ****m, int n_box, int n_page, int n_row) {
  int i, j, k;
  for (i=0; i<n_box; i++) {
    for (j=0; j<n_page; j++) {
      for (k=0; k<n_row; k++) {
        delete m[i][j][k];
      }
      delete m[i][j];
    }
    delete m[i];
  }
  delete m;
}

//...5D array allocation...
double *****matrix5d(int n_cn, int n_box, int n_page, int n_row, int n_col){
  int i, j, k, n, l;
  double *****m;
  m = new double ****[n_cn];
  if (!m) cout << "Insufficient memory in matrix5d()\n";
  for (i=0; i<n_cn; i++) {
    m[i] = new double ***[n_box];
    if (!m[i]) cout << "Insufficient memory in matrix5d()\n";
	  for (j=0; j<n_box; j++) {
      m[i][j] = new double **[n_page];
      if (!m[i][j]) cout << "Insufficient memory in matrix5d()\n";
      for (k=0; k<n_page; k++) {
        m[i][j][k] = new double *[n_row];
				if (!m[i][j][k]) cout << "Insufficient memory in matrix5d()\n";
        for (n=0; n<n_row; n++) {
          m[i][j][k][n] = new double[n_col];
          if (!m[i][j][k][n]) 
              cout << "Insufficient memory in matrix5d()\n";
          for (l=0; l<n_col; l++) {
            m[i][j][k][n][l] = 0.0;
          }
        }//...n
      }//...k
    }//...j
  }//...i
  return m;
}

void free_matrix5d(double *****m, int n_cn, int n_box, int n_page, int n_row){
  int i, j, k, n;
  for (i=0; i<n_cn; i++) {
    for (j=0; j<n_box; j++) {
      for (k=0; k<n_page; k++) {
        for (n=0; n<n_row; n++) {
          delete m[i][j][k][n];
        }
        delete m[i][j][k];
      }
      delete m[i][j];
    }
    delete m[i];
  }
  delete m;
}

