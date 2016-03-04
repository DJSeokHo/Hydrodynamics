//************************************************************************//
//***************************    201150237    ****************************//
//***************************     SHI HAO     ****************************//
//***************************     Seok-Ho     ****************************//
//************************************************************************//

#include < stdio.h >
#include < math.h >
#include < stdlib.h >

#define NP 100
#define SUB 0.5
#define sub SUB + 1
#define thick 0.12
#define pi 3.141592

#define alpha 10

#define arc 2 * PI / NP

void force( void );
void getpot( int NNP, int k, double dblet[ NP * 3 ], double source[ NP * 3 ] );
void matrixsolution( double A[ NP * 3 ][ NP * 3 ], double known[ NP * 3 ], int num_var );
double yt( double , double );

double Xnodp[ NP * 2 ], Ynodp[ NP * 2 ], Xcolp[ NP * 2 ], Ycolp[ NP * 2 ];
double subXnodp[ NP * 2 ], subYnodp[ NP * 2 ], subXcolp[ NP * 2 ], subYcolp[ NP * 2 ];
double zXnodp[ NP * 3 ], zYnodp[ NP * 3 ], zXcolp[ NP * 3 ], zYcolp[ NP * 3 ];

double Plnth[ NP * 2 ], subPlnth[ NP * 2 ], zPlnth[ NP * 3 ];
double phi[ NP * 3 ];
double Nx[ NP * 2 ], Ny[ NP * 2 ], Tx[ NP * 2 ], Ty[ NP * 2 ];
double subNx[ NP * 2 ], subNy[ NP * 2 ], subTx[ NP * 2 ], subTy[ NP * 2 ];
double zNx[ NP * 2 ], zNy[ NP * 2 ], zTx[ NP * 2 ], zTy[ NP * 2 ];

double Vs[ NP * 3 ], Cp[ NP * 3 ];

main()
{
	FILE *fp_geo, *fp_pot, *fp_vel, *fp_pre;
	int i, j, self = 0;
	double xn[ NP * 2 ], yn[ NP * 2 ];
	double subxn[ NP * 2 ], subyn[ NP * 2 ];
	double nceta;
	double Aij[ NP * 3 ][ NP * 3 ], Dblet[ NP * 3 ], Bij[ NP * 3 ][ NP * 3 ], cd[ NP * 3 ][ NP * 3 ], zSigma[ NP * 3 ], RHS[ NP * 3 ];
	double Source[ NP * 3];
	double subSigma[ NP * 3 ], Sigma[ NP * 3 ];
	double indvs[ NP * 3 ];
	double Uxinf, Uyinf, dx1, dx2;
	double kceta, subkceta, KUTTA, subKUTTA, Cl1,Cl2;
	double rxkt, rykt, subrxkt, subrykt;
	for( i = 0; i <= NP * 3; i++ )
	{
		Dblet[ i ] = 0;
		Source[ i ] = 0;
		RHS[ i ] = 0.0;
		Cp[ i ] = 0.0;
	}
	fp_geo = fopen("geo.dat","wt");
	fp_pot = fopen("pot.dat","wt");
	fp_vel = fopen("vel.dat","wt");
	fp_pre = fopen("pre.dat","wt");
	for( i = 1; i <= NP; i++ )
	{
		Xnodp[ i ] = Ynodp[ i ] = 0.0;
		subXnodp[ i ] = subYnodp[ i ] = 0.0;
    }
    for( i = 1; i <= NP + 1; i++ )
	{
		nceta = ( i - 1 ) * arc;
		xn[ i ] = 0.5 * ( 1.0 + cos( nceta ) );
		subxn[ i ] = xn[ i ] + sub;
    }
    Xnodp[ 1 ] = 1.0;
    Ynodp[ 1 ] = 0.0;
	subXnodp[ 1 ] = Xnodp[ 1 ] + sub;
	subYnodp[ 1 ] = Ynodp[ 1 ];
    yn[ 1 ] = 0.0;
    yn[ NP / 2 + 1 ] = 0.0;
   	subyn[ 1 ] = 0.0;
	subyn[ NP / 2 + 1 ] = 0.0;
    for( i = 2; i <= ( NP / 2 ); i++ )
	{
		yn[ i ] = yt( xn[ i ], thick );
		subyn[ i ] = yn[ i ];
		Xnodp[ i ] = xn[ i ];
		Ynodp[ i ] = yn[ i ];
		subXnodp[ i ] = subxn[ i ];
		subYnodp[ i ] = subyn[ i ];
    }
	for(i = ( NP / 2+1 ); i <= NP + 1; i++ )
	{
		Xnodp[ i ] = xn[ i ];
		Ynodp[ i ] = -yn[ NP + 2 - i ];
		subXnodp[ i ] = subxn[ i ];
		subYnodp[ i ] = -subyn[ NP + 2 - i ];
	}
	//**********************************suo you dian zuo biao*******************************
	fprintf( fp_geo, "title =\"thick=%3.2f  Np=%d \"\n", thick, NP );
	fprintf( fp_geo, "zone t=\"Node point\"\n" );
	for( i = 1; i <= NP + 1; i++ )
	{
		fprintf(fp_geo, "%0.5f  %0.5f\n", Xnodp[ i ], Ynodp[ i ] );
	}
	fprintf( fp_geo, "zone t=\"Control point\"\n" );
	for( i = 1; i <= NP; i++ )
	{
		Xcolp[ i ] = ( Xnodp[ i + 1 ] + Xnodp[ i ] )*0.5;
		Ycolp[ i ] = ( Ynodp[ i + 1 ] + Ynodp[ i ] )*0.5;
		fprintf(fp_geo, "%0.5f  %0.5f\n", Xcolp[ i ], Ycolp[ i ] );
		Plnth[i] = sqrt( pow(Xnodp[ i + 1 ] - Xnodp[ i ], 2 ) + pow( Ynodp[ i + 1 ] - Ynodp[ i ], 2 ) );
		Nx[ i ] = ( Ynodp[ i + 1 ] - Ynodp[ i ] ) / Plnth[ i ];
		Ny[ i ] = ( Xnodp[ i ] - Xnodp[ i + 1 ] ) / Plnth[ i ];
		Tx[ i ] = -Ny[ i ]; 
		Ty[ i ] = Nx[ i ];
	}
	//***********************************suo you dian zuo biao*******************************
	//**************************************SUB suo you dian zuo biao************************
	fprintf( fp_geo, "title =\"thick=%3.2f  Np=%d \"\n", thick, NP );
	fprintf( fp_geo, "zone t=\"Node point\"\n" );
	for( i = 1; i <= NP + 1; i++ )
	{
		fprintf( fp_geo, "%0.5f  %0.5f\n", subXnodp[ i ], subYnodp[ i ] );
	}
	fprintf( fp_geo, "zone t=\"Control point\"\n" );
	for( i = 1; i <= NP; i++ )
	{
		subXcolp[ i ] = ( subXnodp[ i + 1 ] + subXnodp[ i ] )*0.5; 
		subYcolp[ i ] = ( subYnodp[ i + 1 ] + subYnodp[ i ] )*0.5;
		fprintf( fp_geo, "%0.5f  %0.5f\n", subXcolp[ i ], subYcolp[ i ] );
		subPlnth[i] = sqrt( pow( subXnodp[ i + 1 ] - subXnodp[ i ], 2 ) + pow( subYnodp[ i + 1 ] - subYnodp[ i ], 2 ) );
		subNx[ i ] = ( subYnodp[ i + 1 ] - subYnodp[ i ] ) / subPlnth[ i ];
		subNy[ i ] = ( subXnodp[ i ] - subXnodp[ i + 1 ] ) / subPlnth[ i ];
		subTx[ i ] = -subNy[ i ]; 
		subTy[ i ] = subNx[ i ];
	}
	//**************************************SUB suo you dian zuo biao************************
	for( i = 1; i <= NP + 1; i++ )
	{
		zXnodp[ i ] = Xnodp[ i ];
		zYnodp[ i ] = Ynodp[ i ];
	}
	for( i = NP + 2, j = 1; i <= ( 2 * NP ) + 3, j <= NP + 1; i++, j++ )
	{
		zXnodp[ i ] = subXnodp[ j ];
		zYnodp[ i ] = subYnodp[ j ];
	}
	for( i = 1; i <= NP; i++ )
	{
		zXcolp[ i ] = Xcolp[ i ];
		zYcolp[ i ] = Ycolp[ i ];
		zPlnth[ i ] = Plnth[ i ];
		zNx[ i ] = Nx[ i ];
		zNy[ i ] = Ny[ i ];
		zTx[ i ] = Tx[ i ];
		zTy[ i ] = Ty[ i ];
	}
	for( i = NP + 1, j = 1; i <= 2 * NP, j <= NP; i++, j++ )
	{
		zXcolp[ i ] = Xcolp[ j ] + sub;
		zYcolp[ i ] = subYcolp[ j ];
		zPlnth[ i ] = subPlnth[ j ];
		zNx[ i ] = subNx[ j ];
		zNy[ i ] = subNy[ j ];
		zTx[ i ] = subTx[ j ];
		zTy[ i ] = subTy[ j ];
	}
	Uxinf = cos( alpha * pi / 180.0 ); /*free stream velocity*/
	Uyinf = sin( alpha * pi / 180.0 );
	for( i = 1; i <= NP; i++ )
	{
		Sigma[ i ] = ( Nx[ i ] * Uxinf + Ny[ i ] * Uyinf );
		subSigma[ i ] = ( subNx[ i ] * Uxinf + subNy[ i ] * Uyinf );
		zSigma[ i ] = Sigma[ i ];
	}
	for( i = NP + 1, j = 1; i <= 2 * NP, j <= NP; i++, j++ )
	{
		zSigma[ i ] = subSigma[ j ];
	}
//*******************************************************setmat**************************************************************
//Page 9 - 11
	for( i = 1; i <= 2 * NP; i++ )
	{
		self = i;
		getpot( 2 * NP, self, Dblet, Source );
		for( j = 1; j <= 2 * NP; j++ )
		{
			cd[ j ][ i ] = 0.0;
			cd[ j ][ j ] = 1.0;
			Aij[ j ][ i ] = cd[ j ][ i ] - Dblet[ j ];
			Bij[ j ][ i ] = Source[ j ];
			RHS[ j ] = RHS[ j ] - zSigma[ i ] * Bij[ j ][ i ];
		}
	}
//*******************************************************setmat**************************************************************
//*******************************************************KUTTA***************************************************************
//Page17
	rxkt =  zXcolp[ 1 ] - zXcolp[ NP ];
	rykt =  zYcolp[ 1 ] - zYcolp[ NP ];

	for( j = 1; j <= 2 * NP; j++ )
	{
		if( j <= NP )
		{
			kceta = atan2( zYcolp[ j ], zXcolp[ j ] - 1.0 );
		}
		else
		{
			kceta = atan2( zYcolp[ j ], 1.0 - zXcolp[ j ] );
		}

		if( kceta < 0.0 )
		{
			KUTTA = ( kceta + 2.0 * pi ) / ( 2.0 * pi );
		}
		else
		{
			KUTTA = kceta / ( 2.0 * pi );
		}
		Aij[ j ][ 1 ] = Aij[ j ][ 1 ] + KUTTA;
		Aij[ j ][ NP ] = Aij[ j ][ NP ] - KUTTA; 
		RHS[ j ] = RHS[ j ] - ( Uxinf * rxkt + Uyinf * rykt ) * KUTTA;
	}
	subrxkt =  zXcolp[ NP + 1 ] - zXcolp[ 2 * NP ];
	subrykt =  zYcolp[ NP + 1 ] - zYcolp[ 2 * NP ];
	for( j = 1; j <= 2 * NP; j++ )
	{
		subkceta = atan2( zYcolp[ j ], zXcolp[ j ] - ( sub + 1 ) );
		if( subkceta < 0.0 )
		{
			subKUTTA = ( subkceta + 2.0 * pi ) / ( 2.0 * pi );
		}
		else
		{
			subKUTTA = subkceta / ( 2.0 * pi );
		}
		Aij[ j ][ NP + 1] = Aij[ j ][ NP + 1] + subKUTTA;
		Aij[ j ][ 2 * NP ] = Aij[ j ][ 2 * NP ] - subKUTTA;
		RHS[ j ] = RHS[ j ] - ( Uxinf * subrxkt + Uyinf * subrykt ) * subKUTTA;
	}
	matrixsolution( Aij, RHS, 2 * NP );
	phi[ 0 ] = phi[ 1 ];
	phi[ 2 * NP + 1 ] = phi[ 2 * NP ];
	Cl1 = 2.0 * ( phi[ 1 ] - phi[ NP ] ); 
	Cl2 =2.0 * ( phi[ NP+1 ] - phi[ 2*NP ] );
	 fprintf(fp_pot,"title=\"t=%2.2f %dPANEL alpha= %d \"\n",thick,NP,alpha);
	 fprintf(fp_vel,"title=\"t=%2.2f %dPANEL alpha= %d \"\n",thick,NP,alpha);
	 fprintf(fp_pre,"title=\"t=%2.2f %dPANEL alpha= %d \"\n",thick,NP,alpha);
	 fprintf(fp_pot,"variables=\"x\",\"y\"\n");
	 fprintf(fp_vel,"variables=\"x\",\"y\"\n");
	 fprintf(fp_pre,"variables=\"x\",\"y\"\n");
	 fprintf(fp_pot,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl1);
	 fprintf(fp_vel,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl1);
	 fprintf(fp_pre,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl1);
	for( i = 1; i <= NP; i++ )
	{
		if( i >= 1 && i <= 3 )
		{
		  dx1 = 0.5 * ( zPlnth[ i ] + zPlnth[ i + 1 ] );
		  dx2 = 0.5 * ( zPlnth[ i + 1 ] + zPlnth[ i + 2 ] );
		  indvs[ i ] = -( dx1 * dx1 * phi[ i + 2 ] -
			 ( dx1 + dx2 ) * ( dx1 + dx2 ) * phi[ i + 1 ] +
			 ( 2. * dx1 * dx2 + dx2 * dx2 ) * phi[ i ] ) /
			 ( dx1 * dx2 * ( dx1 + dx2 ) );
		}
		else if( i >= NP - 2 && i <= NP )
		{
		  dx1 = 0.5 * ( zPlnth[ i ] + zPlnth[ i - 1 ] );
		  dx2 = 0.5 * ( zPlnth[ i - 1 ] + zPlnth[ i - 2 ] );
		  indvs[ i ] = ( dx1 * dx1 * phi[ i - 2 ] -
			( dx1 + dx2 ) * ( dx1 + dx2 ) * phi[ i - 1 ] + 
			( 2. * dx1 * dx2 + dx2 * dx2 ) * phi[ i ] ) /
			( dx1 * dx2 * ( dx1 + dx2 ) );
		}
		else
		{
		  dx1 = 0.5 * ( zPlnth[ i ] + zPlnth[ i - 1 ] );
		  dx2 = 0.5 * ( zPlnth[ i ] + zPlnth[ i + 1 ] );
		  indvs[ i ] = ( dx1 * dx1 * ( phi[ i + 1 ] - phi[ i ] ) +
			dx2 * dx2 * ( phi[ i ] - phi[ i - 1 ] ) ) /
			( dx1 * dx2 * ( dx1 + dx2 ) );
		}
		Vs[ i ] = -( Uxinf * zTx[ i ] + Uyinf * zTy[ i ] + indvs[ i ] );
		Vs[ 1 ] = Vs[ NP ] = 0.0;
		Cp[ i ] = 1.0 - ( Vs[ i ] * Vs[ i ] );
		fprintf(fp_pot,"%3.6f  %3.6f  \n",zXcolp[i],phi[i]);
		fprintf(fp_vel,"%3.6f  %3.6f  \n",zXcolp[i],fabs(Vs[i]));
		fprintf(fp_pre,"%3.6f  %3.6f  \n",zXcolp[i],-Cp[i]);
	}
	fprintf(fp_pot,"title=\"t=%2.2f %dPANEL alpha= %d \"\n",thick,NP,alpha);
	fprintf(fp_vel,"title=\"t=%2.2f %dPANEL alpha= %d \"\n",thick,NP,alpha);
	fprintf(fp_pre,"title=\"t=%2.2f %dPANEL alpha= %d \"\n",thick,NP,alpha);
	fprintf(fp_pot,"variables=\"x\",\"y\"\n");
	fprintf(fp_vel,"variables=\"x\",\"y\"\n");
	fprintf(fp_pre,"variables=\"x\",\"y\"\n");
	fprintf(fp_pot,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl2);
	fprintf(fp_vel,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl2);
	fprintf(fp_pre,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl2);
	for( i = NP + 1; i <= 2 * NP; i++ )
	{
		if( i >= NP + 1 && i <= NP + 3 )
		{
		  dx1 = 0.5 * ( zPlnth[ i ] + zPlnth[ i + 1 ] );
		  dx2 = 0.5 * ( zPlnth[ i + 1 ] + zPlnth[ i + 2 ] );
		  indvs[ i ] = -( dx1 * dx1 * phi[ i + 2 ] -
			 ( dx1 + dx2 ) * ( dx1 + dx2 ) * phi[ i + 1 ] +
			 ( 2. * dx1 * dx2 + dx2 * dx2 ) * phi[ i ] ) /
			 ( dx1 * dx2 * ( dx1 + dx2 ) );
		}
		else if( i >= 2 * NP - 2 && i <= 2 * NP )
		{
		  dx1 = 0.5 * ( zPlnth[ i ] + zPlnth[ i - 1 ] );
		  dx2 = 0.5 * ( zPlnth[ i - 1 ] + zPlnth[ i - 2 ] );
		  indvs[ i ] = ( dx1 * dx1 * phi[ i - 2 ] -
			( dx1 + dx2 ) * ( dx1 + dx2 ) * phi[ i - 1 ] + 
			( 2. * dx1 * dx2 + dx2 * dx2 ) * phi[ i ] ) /
			( dx1 * dx2 * ( dx1 + dx2 ) );
		}
		else
		{
		  dx1 = 0.5 * ( zPlnth[ i ] + zPlnth[ i - 1 ] );
		  dx2 = 0.5 * ( zPlnth[ i ] + zPlnth[ i + 1 ] );
		  indvs[ i ] = ( dx1 * dx1 * ( phi[ i + 1 ] - phi[ i ] ) +
			dx2 * dx2 * ( phi[ i ] - phi[ i - 1 ] ) ) /
			( dx1 * dx2 * ( dx1 + dx2 ) );
		}
		Vs[ i ] = -( Uxinf * zTx[ i ] + Uyinf * zTy[ i ] + indvs[ i ] );
		Vs[ NP + 1 ] = Vs[ 2 * NP ] = 0.0;
		Cp[ i ] = 1.0 - ( Vs[ i ] * Vs[ i ] );
		fprintf(fp_pot,"%3.6f  %3.6f  \n",zXcolp[i],phi[i]);
		fprintf(fp_vel,"%3.6f  %3.6f  \n",zXcolp[i],fabs(Vs[i]));
		fprintf(fp_pre,"%3.6f  %3.6f  \n",zXcolp[i],-Cp[i]);
	}
	force();
	return 0;
}
//*******************************************************KUTTA***************************************************************
void force( void )
{
	int i;
	double Fy, Fx, Fxd, Fd, cf, subFy, subFx, subFxd, subFd, subcf;
	double CL, CD, CLf, CDf, subCL, subCD, subCLf, subCDf;
	Fy = Fx = Fd = Fxd = subFy = subFx = subFd = subFxd = CL = CLf = CD = CDf = subCL = subCLf = subCD = subCDf = 0.0;
	cf = subcf = 0.0035;
	for( i = 1; i <= NP; i++ )
	{
		Fy += -Cp[ i ] * zNy[ i ] * zPlnth[ i ];
		Fx += -Cp[ i ] * zNx[ i ] * zPlnth[ i ];
		Fd += cf * Vs[ i ] * Vs[ i ] * zPlnth[ i ];
		Fxd = Fx + Fd;
		CL = Fy * cos( alpha * pi / 180.0 ) - Fx * sin( alpha * pi / 180.0 );
		CD = Fx * cos( alpha * pi / 180.0 ) + Fy * sin( alpha * pi / 180.0 );
 		CLf = Fy * cos( alpha * pi / 180.0 ) - Fxd * sin( alpha * pi / 180.0 );
		CDf = Fxd * cos( alpha * pi / 180.0 ) + Fy * sin( alpha * pi / 180.0 );
	}
		printf("\n=========================Forw===========================\n");
		printf(" alpha     CL         CD          CLf         CDf\n");
		printf("----------------------------------------------------------\n");
		printf( "   %d  ", alpha );
		printf( "   %0.5lf  ", CL );
		printf( "   %0.5lf  ", CD );
		printf( "   %0.5lf  ", CLf );
		printf( "   %0.5lf  \n", CDf );
		printf("\n----------------------------------------------------------\n");
		printf("\n==========================================================\n");

	for( i = NP + 1; i <= 2 * NP; i++ )
	{
		subFy += -Cp[ i ] * zNy[ i ] * zPlnth[ i ];
		subFx += -Cp[ i ] * zNx[ i ] * zPlnth[ i ];
		subFd += subcf * Vs[ i ] * Vs[ i ] * zPlnth[ i ];
		subFxd = subFx + subFd;
  		subCL = subFy * cos( alpha * pi / 180.0 ) - subFx * sin( alpha * pi / 180.0 );
		subCD = subFx * cos( alpha * pi / 180.0 ) + subFy * sin( alpha * pi / 180.0 );
 		subCLf = subFy * cos( alpha * pi / 180.0 ) - subFxd * sin( alpha * pi / 180.0 );
		subCDf = subFxd * cos( alpha * pi / 180.0 ) + subFy * sin( alpha * pi / 180.0 );
	}
		printf("\n=========================Aft==============================\n");
		printf(" alpha     CL         CD          CLf         CDf\n");
		printf("----------------------------------------------------------\n");
		printf( "   %d  ", alpha );
		printf( "   %0.5lf  ", subCL );
		printf( "   %0.5lf  ", subCD );
		printf( "   %0.5lf  ", subCLf );
		printf( "   %0.5lf  \n", subCDf );
		printf("\n----------------------------------------------------------\n");
		printf("\n==========================================================\n");
}

void getpot( int NNP, int k, double dblet[ NP * 3 ], double source[ NP * 3 ] ) //Page 10
{
	  int j;
	  double tx, ty, xn1, yn1, xn2, yn2, xc, yc, pl;
	  double l1, l2, twr1, twr2, h, theta;
	  if( k <= NNP / 2 )
	  {
		  xn1 = zXnodp[ k ]; // node point
		  yn1 = zYnodp[ k ]; // node point
		  xn2 = zXnodp[ k + 1 ]; // node point
		  yn2 = zYnodp[ k + 1 ]; // node point
		  pl = zPlnth[ k ]; //panel length
	  }
	  else
	  {
		  xn1 = zXnodp[ k + 1]; // node point
		  yn1 = zYnodp[ k + 1 ]; // node point
		  xn2 = zXnodp[ k + 2 ]; // node point
		  yn2 = zYnodp[ k + 2 ]; // node point
		  pl = zPlnth[ k ]; //panel length
	  }
		  tx = ( xn2 - xn1 ) / pl;
		  ty = ( yn2 - yn1 ) / pl;
	  for( j = 1; j <= NNP; j++ )
	  {
		xc = zXcolp[ j ]; // control point
		yc = zYcolp[ j ]; // control point
		l1 = ( xc - xn1 ) * tx + ( yc - yn1 ) * ty;
		l2 = -( ( xc - xn2 ) * tx + ( yc - yn2 ) * ty );
		h = -( ( xc - xn1 ) * ty - ( yc - yn1 ) * tx );
		twr1 = ( xc - xn1 ) * ( xc - xn1 ) + ( yc - yn1 ) * ( yc - yn1 ); /*squre(r1)*/
		twr2 = ( xc - xn2 ) * ( xc - xn2 ) + ( yc - yn2 ) * ( yc - yn2 ); /*squre(r2)*/
		theta = atan2( fabs( h ) * pl, twr1 - pl * l1 ); //Dblet
		if( j == k )
		{
		  dblet[ j ] = 0.5;
		  source[ j ] = 0.0;
		}
		else
		  dblet[ j ] = -h / fabs( h ) * theta / ( 2.0 * pi );
		  source[ j ] = 0.5 * ( ( pl - l1 ) * log( twr2 ) + l1 * log( twr1 ) - 2.0 * pl + 2.0 * fabs( h ) * theta ) / ( 2.0 * pi );
	  }
}
void matrixsolution( double A[ NP * 3 ][ NP * 3 ], double known[ NP * 3 ], int num_var )
{
	int num_eq;
	int k, l;
	int row, column, pivot_column, max_index;
	double max_value, ftemp1, ftemp2, pivot_value;
	double coeff[ NP * 3 ][ NP * 3 ], inv_coeff[ NP * 3 ][ NP * 3 ], var[ NP * 3 ];
	num_eq = num_var;
 	for( k = 1; k <= num_var; k++ )
	{
		for( l = 1; l <= num_var; l++ )
		{
			coeff[ k ][ l ] = A[ k ][ l ];
		}
	}
	for( row = 1; row <= num_eq; row++ )
		for( column = 1; column <= num_eq; column++ ){
			if( row == column )
				inv_coeff[row][column] = 1;
			else
				inv_coeff[row][column] = 0;
		}
		for( pivot_column = 1; pivot_column <= num_eq; pivot_column++ ){
			max_index = coeff[0][column];
			max_value = 0;
			for( row = pivot_column; row <= num_eq; row++ )
				if( coeff[row][pivot_column]*coeff[row][pivot_column] > max_value*max_value ){
					max_index = row;
					max_value = coeff[row][pivot_column];
				}
				if( pivot_column != max_index )
					for( column = 1; column <= num_eq; column++ ){
						ftemp1 = coeff[pivot_column][column];
						ftemp2 = inv_coeff[pivot_column][column];
						coeff[pivot_column][column] = coeff[max_index][column];	
						inv_coeff[pivot_column][column] = inv_coeff[max_index][column];
						coeff[max_index][column] = ftemp1;
						inv_coeff[max_index][column] = ftemp2;
					}
					pivot_value = coeff[pivot_column][pivot_column];
					for(column = 1; column <= num_eq; column++ ){
						coeff[pivot_column][column] /= pivot_value;
						inv_coeff[pivot_column][column] /= pivot_value;
						}
						for( row = 1; row <= num_eq; row++ )
							if( row != pivot_column ){
								ftemp1 = coeff[row][pivot_column];
								for( column = 1; column <= num_eq; column++ ){
									coeff[row][column] -= ftemp1*coeff[pivot_column][column];
									inv_coeff[row][column] -= ftemp1*inv_coeff[pivot_column][column];
								}
							}
	}
	for( k = 1; k <= num_var; k++ )
	{
		var[ k ] = 0.0;
		for( l = 1; l <= num_eq; l++ )
			var[ k ] += inv_coeff[ k ][ l ] * known[ l ];
	}
	for( k = 1; k <= num_var; k++ )
	{
		phi[ k ] = var[ k ];
	}
}

double yt(double x, double t)
{
	double y=0.0;
	if(x==1.0)
		return y=0.0;
	if(x==0.0)
		return y=0.0;
	else 
	{
		y=(t/0.2)*(0.2969*pow(x,0.5)-0.126*x-0.3537*x*x
		   +0.2843*pow(x,3.)-0.1015*pow(x,4.));
    return y;
	}
}
