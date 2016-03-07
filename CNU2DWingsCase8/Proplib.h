/*----------------proplib.h-------------------------*/
// #define max(i,j)  ( ((i)>(j)) ? (i) : (j)) 
using namespace std;

#define fmax(x,y) ( ((x)>(y)) ? (x) : (y) )
#define fmin(x,y) ( ((x)<(y)) ? (x) : (y) )

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

void drivdk(int,int,double*,double*,double*,double*,double*);
void evaldk(int,int,double*,double*,double*,double*);
double fillin(double,double*,double*,int);
/*
void intedk(int nin, double *Xin, double xl, double xu,
            double &ydx, double &xydx, double &xxydx, double *A);
*/
void intedk(int,double*,double,double,double*,double*,double*,double*);
void simq(double*,double*,int,int*);
void uglydk(int,int,int,double*,double*,double,double,double*);
void iterlu(double**,double**,int,int*,double*,double*,int,double);
void linsol(double**,double**,int,int*,double*,double*,
        int,int,int,double);
void lubksb(double**,int,int*,double*);
void ludcmp(double**,double**,int,int*,double*);
void vorseg(double, double, double, double, double, double,
            double, double, double, double*, double*, double*,
            double*, double*, double*, int);
void prcros(double[], double[], double[]);
double prdot(double[], double[]);
void prcrosK(double *a, double *b, double *c);
double prdotK(double *a, double *b);
double potder(double xp, double xin[], double yin[]);
/*--Exponential spline--*/
void spl(int npt, double Xin[], double Yin[], double *slpl, double *slpr,
         double relax, int itmax, int ispl, double R[], double P[], int *iter);
void espl(int npt, double Xin[], double Yin[], int ispl, double F2[],
          double P[], int nout, double Xout[], double Yout[]);
void dspl(int npt, double Xin[], double Yin[], int ispl, double F2[],
          double P[],  int nout, double Xout[], double Y1out[], double Y2out[]);
void intspl(int npt, double Xin[], double Yin[], int ispl,
            double F2[], double P[], double xa, double xb, double *yint);
void tridag(int il, int iu, double A[], double B[], double C[], double R[]);

void splK(int ispl, double relax, double slpl, double slpr, 
          int npt, double *Xin, double *Yin, double *R, double *P);
void esplK(int ispl, int npt, double *Xin, double *Yin, 
           int nout, double *Xout, double *Yout, double *F2, double *P);
void dsplK(int ispl, int npt, double *Xin, double *Yin,
           int nout, double *Xout, double *Y1out, double *Y2out, double *F2, double *P);
void intsplK(int ispl, int npt, double *Xin, double *Yin, 
            double xa, double xb, double &yint, double *F2, double *P);

//--Numerical Recipes Utilities Header--
double powi(double x, int n);
double powr(double x, double apower); 
double ran3(int *);
double *vector(int,int);
double **matrix(int,int,int,int);
int *ivector(int,int);
void free_vector(double *,int,int);
void free_ivector(int *,int,int);
void free_matrix(double **,int,int,int,int);
void nrerror(char error_text[]);
                                  
//...Memory allocation routines........
int *ivec(int n_row);
double *vec(int n_row);
double **matrix2d(int n_row, int n_col);
void free_matrix2d(double **m, int n_row);
int **imatrix2d(int n_row, int n_col);
void free_imatrix2d(int **m, int n_row);
double ***matrix3d(int n_page, int n_row, int n_col);
void free_matrix3d(double ***m, int n_page, int n_row);
double ****matrix4d(int n_box, int n_page, int n_row, int n_col);
void free_matrix4d(double ****m, int n_box, int n_page, int n_row);
double *****matrix5d(int n_cn, int n_box, int n_page, int n_row, int n_col);
void free_matrix5d(double *****m, int n_cn, int n_box, int n_page,int n_row);

void displayright(void)
{
	cout << "************************************************" << endl;
	cout << "*                CNU2D  case8                  *" << endl;
	cout << "*                                              *" << endl;
	cout << "*                   SeokHo                     *" << endl;
	cout << "*                                              *" << endl;
	cout << "************************************************" << endl;
}



