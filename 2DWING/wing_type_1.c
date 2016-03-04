#include < stdio.h >
#include < math.h >
#include < stdlib.h >

#define NP 100
#define SUB 0.5
#define sub SUB + 1
#define thick 0.12
#define PI 3.141592
#define pi 3.141592
#define Alpha 0 

#define alpha 0 // yang jiao
#define hudu -alpha * PI / 180
#define arc 2 * PI / NP
#define rx 1
#define ry 0 

#define betax 0.850   
#define betay 0  
 
#define beta 0 
#define betahudu -beta * PI / 180

void force(void);
void getpot( int NNP, int k, double dblet[ 200 ], double source[ 200 ] );
void matrixsolution( double A[ 200 ][ 200 ], double known[ 200 ], int num_var );

double realsin( double );
double realcos( double );
double yc(double, double, double);
double yt(double , double );
double RotaryAngleX( double, double );
double RotaryAngleY( double, double );

//**************************   hudu   ***************************************************************************
double BetaRotaryAngleX( double, double );
double BetaRotaryAngleY( double, double );


double Xnodp[ 200 ], Ynodp[ 200 ], Xcolp[ 200 ], Ycolp[ 200 ];
double subXnodp[ 200 ], subYnodp[ 200 ], subXcolp[ 200 ], subYcolp[ 200 ];
double Plnth[ 200 ];
double phi[ 200 ];
double Nx[ 200 ], Ny[ 200 ], Tx[ 200 ], Ty[ 200 ];
double Vs[ 200 ], Cp[ 200 ];

main()
{
	FILE *fp_geo, *fp_pot, *fp_vel, *fp_pre;

	int i, j, self = 0;
	
	double xn[ 200 ], yn[ 200 ];
	double nceta;
	
	double RAXUP[ 200 ], RAYUP[ 200 ], RAXDOWN[ 200 ], RAYDOWN[ 200 ], RAXcolp[ 200 ], RAYcolp[ 200 ];
	double BRAXUP[ 200 ], BRAYUP[ 200 ], BRAXDOWN[ 200 ], BRAYDOWN[ 200 ];
	double BRAXcolpUP[ 200 ], BRAYcolpUP[ 200 ], BRAXcolpDOWN[ 200 ], BRAYcolpDOWN[ 200 ];
	double Dblet[ 200 ], Source[ 200 ];
	double Aij[ 200 ][ 200 ], Bij[ 200 ][ 200 ], cd[ 200 ][ 200 ], Sigma[ 200 ], RHS[ 200 ];
	double indvs[ 200 ];
	double Uxinf, Uyinf, dx1, dx2;
	double kceta, KUTTA, Cl;
	double rxkt, rykt;

	for( i = 0; i <= 200; i++ )
	{
		Dblet[ i ] = 0;
		Source[ i ] = 0;
	}

	fp_geo = fopen("geo.dat","wt");
	fp_pot = fopen("pot.dat","wt");
	fp_vel = fopen("vel.dat","wt");
	fp_pre = fopen("pre.dat","wt");

	for( i = 1; i <= NP; i++ )
	{
		Xnodp[ i ] = Ynodp[ i ] = 0.0;
    }
    
    for( i = 1; i <= NP + 1; i++ )
	{
		nceta = ( i - 1 ) * arc;
		xn[ i ] = 0.5 * ( 1.0 + cos( nceta ) );
    }

    Xnodp[ 1 ] = 1.0;
    Ynodp[ 1 ] = 0.0;
    yn[ 1 ] = 0.0;
    yn[ NP / 2 + 1 ] = 0.0;
    
    for( i = 1; i <= ( NP / 2 ); i++ )
	{
		yn[ i ] = yt( xn[ i ], thick );
		Xnodp[ i ] = xn[ i ];
		Ynodp[ i ] = yn[ i ];
	
		RAXUP[ i ] = RotaryAngleX( Xnodp[ i ], Ynodp[ i ] );
		RAYUP[ i ] = RotaryAngleY( Xnodp[ i ], Ynodp[ i ] );
    }
	
	//*******************************up body*****************


	//printf("\n***************************************\n");

    for(i = ( NP / 2+1 ); i <= NP + 1; i++ )
	{
		Xnodp[ i ] = xn[ i ];
		Ynodp[ i ] = -yn[ NP + 2 - i ];
	
		RAXDOWN[ i ] = RotaryAngleX( Xnodp[ i ], Ynodp[ i ] );
		RAYDOWN[ i ] = RotaryAngleY( Xnodp[ i ], Ynodp[ i ] );
	}
	//*******************************down body*******************
	
	//****************************************midpoint**************************************
	for( i = 1; i <= NP; i++ )
	{
		Xcolp[ i ] = ( Xnodp[ i + 1 ] + Xnodp[ i ] )*0.5; //轟嵐실櫓듐
		Ycolp[ i ] = ( Ynodp[ i + 1 ] + Ynodp[ i ] )*0.5;
		Plnth[i] = sqrt( pow(Xnodp[ i + 1 ] - Xnodp[ i ], 2 ) + pow( Ynodp[ i + 1 ] - Ynodp[ i ], 2 ) );

		Nx[ i ] = ( Ynodp[ i + 1 ] - Ynodp[ i ] ) / Plnth[ i ];
		Ny[ i ] = ( Xnodp[ i ] - Xnodp[ i + 1 ] ) / Plnth[ i ];
		Tx[ i ] = -Ny[ i ]; /*(Xnodp[i+1]-Xnodp[i])/Plnth[i];*/
		Ty[ i ] = Nx[ i ]; /*(Ynodp[i+1]-Ynodp[i])/Plnth[i];*/

		RAXcolp[ i ] = RotaryAngleX( Xcolp[ i ], Ycolp[ i ] );
		RAYcolp[ i ] = RotaryAngleY( Xcolp[ i ], Ycolp[ i ] );
	}

	for( j = 1; j <= 12; j++ )
	{
		BRAXcolpUP[ j ] = BetaRotaryAngleX( RAXcolp[ j ], RAYcolp[ j ] );
		BRAYcolpUP[ j ] = BetaRotaryAngleY( RAXcolp[ j ], RAYcolp[ j ] );
		Xcolp[ j ] = BRAXcolpUP[ j ];
		Ycolp[ j ] = BRAYcolpUP[ j ];
	}

	for( i = 13; i <= 50; i++ )
	{
		Xcolp[ i ] = RAXcolp[ i ];
		Ycolp[ i ] = RAYcolp[ i ];
	}
	for( i = 51; i <= 88; i++ )
	{
		Xcolp[ i ] = RAXcolp[ i ];
		Ycolp[ i ] = RAYcolp[ i ];
	}

	for( j = 89; j <= 100; j++ )
	{
		BRAXcolpDOWN[ j ] = BetaRotaryAngleX( RAXcolp[ j ], RAYcolp[ j ] );
		BRAYcolpDOWN[ j ] = BetaRotaryAngleY( RAXcolp[ j ], RAYcolp[ j ] );
		
		Xcolp[ j ] = BRAXcolpDOWN[ j ];
		Ycolp[ j ] = BRAYcolpDOWN[ j ];
	}

	for( j = 1; j <= 13; j++ ) 
	{
		BRAXUP[ j ] = BetaRotaryAngleX( RAXUP[ j ], RAYUP[ j ] );
		BRAYUP[ j ] = BetaRotaryAngleY( RAXUP[ j ], RAYUP[ j ] );
		Xnodp[ j ] = BRAXUP[ j ];
		Ynodp[ j ] = BRAYUP[ j ];
	}

	
	for(i = 14;i <= NP / 2;i++)
	{
		Xnodp[ i ] = RAXUP[ i ];
		Ynodp[ i ] = RAYUP[ i ];
	}

	for(i=51;i<=88;i++)
	{
		Xnodp[ i ] = RAXDOWN[ i ];
		Ynodp[ i ] = RAYDOWN[ i ];
	}

//	printf("\n****************BeTA***********************\n");

	for(j=89;j<=101;j++)
	{
		BRAXDOWN[ j ] = BetaRotaryAngleX( RAXDOWN[ j ], RAYDOWN[ j ] );
		BRAYDOWN[ j ] = BetaRotaryAngleY( RAXDOWN[ j ], RAYDOWN[ j ] );

		Xnodp[ j ] = BRAXDOWN[ j ];
		Ynodp[ j ] = BRAYDOWN[ j ];
	//****************************************nodpoint**************************************
	}

//**************************************suo you dian zuo biao*******************************

	fprintf( fp_geo, "title =\"thick=%3.2f  Np=%d \"\n", thick, NP );
	fprintf( fp_geo, "zone t=\"Node point\"\n" );
	for( i = 1; i <= 101; i++ )
	{
		fprintf(fp_geo, "%0.5f  %0.5f\n", Xnodp[ i ], Ynodp[ i ] );
	}

	fprintf( fp_geo, "zone t=\"Control point\"\n" );
	for( i = 1; i <= 100; i++ )
	{
		fprintf(fp_geo, "%0.5f  %0.5f\n", Xcolp[ i ], Ycolp[ i ] );
	}
//**************************************suo you dian zuo biao*******************************

//**************************************SUB suo you dian zuo biao***************************

	for( i = 1; i <= 101; i++ )
	{
		subXnodp[ i ] = Xnodp[ i ] + sub;
		subYnodp[ i ] = Ynodp[ i ];
	}
	for( i = 1; i <= 100; i++ )
	{
		subXcolp[ i ] = Xcolp[ i ] + sub;
		subYcolp[ i ] = Ycolp[ i ];
	}

	fprintf( fp_geo, "title =\"thick=%3.2f  Np=%d \"\n", thick, NP );
	fprintf( fp_geo, "zone t=\"Node point\"\n" );
	for( i = 1; i <= 101; i++ )
	{
		fprintf( fp_geo, "%0.5f  %0.5f\n", subXnodp[ i ], subYnodp[ i ] );
	}

	fprintf( fp_geo, "zone t=\"Control point\"\n" );
	for( i = 1; i <= 100; i++ )
	{
		fprintf( fp_geo, "%0.5f  %0.5f\n", subXcolp[ i ], subYcolp[ i ] );
	}
//**************************************SUB suo you dian zuo biao***************************

	Uxinf = cos( Alpha * pi / 180.0 ); /*free stream velocity*/
	Uyinf = sin( Alpha * pi / 180.0 );

	for( i = 1; i <= NP; i++ )
	{
		RHS[ i ] = 0.0;
		Sigma[ i ] = ( Nx[ i ] * Uxinf + Ny[ i ] * Uyinf );
	}

	Plnth[ 0 ] = Plnth[ 1 ];
	Plnth[ NP + 1 ] = Plnth[ NP ];
//***************************************************************************************************************************
//*******************************************************setmat**************************************************************
//Page 9 - 11
	for( i = 1; i <= NP; i++ )
	{
		self = i;
		getpot( NP, self, Dblet, Source );

		for( j = 1; j <= NP; j++ )
		{
			cd[ j ][ i ] = 0.0;
			cd[ j ][ j ] = 1.0;
		
			Aij[ j ][ i ] = cd[ j ][ i ] - Dblet[ j ];
			Bij[ j ][ i ] = Source[ j ];
			
			RHS[ j ] = RHS[ j ] - Sigma[ i ] * Bij[ j ][ i ];
		}
	}
	
//*******************************************************setmat**************************************************************

//*******************************************************KUTTA***************************************************************
//Page17
	rxkt =  Xcolp[ 1 ] - Xcolp[ NP ];
	rykt =  Ycolp[ 1 ] - Ycolp[ NP ];

	for( j = 1; j <= NP; j++ )
	{
		//***********************************ceta hou mian*****************************
		kceta = atan2( Ycolp[ j ], Xcolp[ j ] - 1.0 );
		if( kceta < 0.0)
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
	matrixsolution( Aij, RHS, NP );
	

	phi[ 0 ] = phi[ 1 ];
	phi[ NP + 1 ] = phi[ NP ];

	
	
	Cl = 2.0 * ( phi[ 1 ] - phi[ NP ] ); 

	for( i = 1; i <= NP; i++ )
	{
		if( i >= 1 && i <= 3 )
		{
		  dx1 = 0.5 * ( Plnth[ i ] + Plnth[ i + 1 ] );
		  dx2 = 0.5 * ( Plnth[ i + 1 ] + Plnth[ i + 2 ] );
		  indvs[ i ] = -( dx1 * dx1 * phi[ i + 2 ] -
			 ( dx1 + dx2 ) * ( dx1 + dx2 ) * phi[ i + 1 ] +
			 ( 2. * dx1 * dx2 + dx2 * dx2 ) * phi[ i ] ) /
			 ( dx1 * dx2 * ( dx1 + dx2 ) );
		}
		else if( i >= NP - 2 && i <= NP )
		{
		  dx1 = 0.5 * ( Plnth[ i ] + Plnth[ i - 1 ] );
		  dx2 = 0.5 * ( Plnth[ i - 1 ] + Plnth[ i - 2 ] );
		  indvs[ i ] = ( dx1 * dx1 * phi[ i - 2 ] -
			( dx1 + dx2 ) * ( dx1 + dx2 ) * phi[ i - 1 ] + 
			( 2. * dx1 * dx2 + dx2 * dx2 ) * phi[ i ] ) /
			( dx1 * dx2 * ( dx1 + dx2 ) );
		}
		else
		{
		  dx1 = 0.5 * ( Plnth[ i ] + Plnth[ i - 1 ] );
		  dx2 = 0.5 * ( Plnth[ i ] + Plnth[ i + 1 ] );
		  indvs[ i ] = ( dx1 * dx1 * ( phi[ i + 1 ] - phi[ i ] ) +
			dx2 * dx2 * ( phi[ i ] - phi[ i - 1 ] ) ) /
			( dx1 * dx2 * ( dx1 + dx2 ) );
		}
		Vs[ i ] = -( Uxinf * Tx[ i ] + Uyinf * Ty[ i ] + indvs[ i ] );
		Vs[ 1 ] = Vs[ NP ] = 0.0;
		Cp[ i ] = 1.0 - ( Vs[ i ] * Vs[ i ] );
	}
	
	  fprintf(fp_pot,"title=\"t=%2.2f %dPANEL alpha= %d \"\n",thick,NP,alpha);
	  fprintf(fp_vel,"title=\"t=%2.2f %dPANEL alpha= %d \"\n",thick,NP,alpha);
	  fprintf(fp_pre,"title=\"t=%2.2f %dPANEL alpha= %d \"\n",thick,NP,alpha);
	  fprintf(fp_pot,"variables=\"x\",\"y\"\n");
	  fprintf(fp_vel,"variables=\"x\",\"y\"\n");
	  fprintf(fp_pre,"variables=\"x\",\"y\"\n");
	  fprintf(fp_pot,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl);
	  fprintf(fp_vel,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl);
	  fprintf(fp_pre,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl);
		
	  for(i=1;i<=NP;i++)
	  {
		fprintf(fp_pot,"%3.6f  %3.6f  \n",Xcolp[i],phi[i]);
		fprintf(fp_vel,"%3.6f  %3.6f  \n",Xcolp[i],fabs(Vs[i]));
		fprintf(fp_pre,"%3.6f  %3.6f  \n",Xcolp[i],-Cp[i]);
	  }
	
	  force();
	return 0;
}
//*******************************************************KUTTA***************************************************************

void force( void )
{
	int i;
	double Fy, Fx, Fxd, Fd, cf;
	double CL, CD, CLf, CDf;

	Fy = Fx = Fd = Fxd = 0.0;
	CL = CLf = CD = CDf = 0.0;
	cf = 0.0035;

	for( i = 1; i <= NP; i++ )
	{
		Fy += -Cp[ i ] * Ny[ i ] * Plnth[ i ];
		Fx += -Cp[ i ] * Nx[ i ] * Plnth[ i ];
		Fd += cf * Vs[ i ] * Vs[ i ] * Plnth[ i ];
	}
    Fxd = Fx + Fd;
  
	CL = Fy * cos( Alpha * pi / 180.0 ) - Fx * sin( Alpha * pi / 180.0 );
	CD = Fx * cos( Alpha * pi / 180.0 ) + Fy * sin( Alpha * pi / 180.0 );
 
	CLf = Fy * cos( Alpha * pi / 180.0 ) - Fxd * sin( Alpha * pi / 180.0 );
	CDf = Fxd * cos( Alpha * pi / 180.0 ) + Fy * sin( Alpha * pi / 180.0 );

	printf("\n==========================================================\n");
	printf(" alpha     CL         CD          CLf         CDf\n");
	printf("----------------------------------------------------------\n");
	printf( "   %d  ", Alpha );
	printf( "   %0.5lf  ", CL );
	printf( "   %0.5lf  ", CD );
	printf( "   %0.5lf  ", CLf );
	printf( "   %0.5lf  \n", CDf );
	printf("\n----------------------------------------------------------\n");
	printf("\n==========================================================\n");
	
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

double yc( double x, double m, double p )
{
  double y;
  if( x <= p ) 
	  return y = m * ( 2 * p * x - x * x ) / ( p * p );
  else 
	  return y = m * ( ( 1 - 2 * p ) + 2 * p * x - x * x ) / ( 1 - 2 * p + p * p );
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
double RotaryAngleX( double RAXX, double RAXY )
{
	double xr;
	xr = ( RAXX - rx ) * cos( hudu ) - ( RAXY - ry ) * sin( hudu ) + rx;
	return xr;
}

double RotaryAngleY( double RAYX, double RAYY )
{
	double yr;
	yr = ( RAYX - rx ) * sin( hudu ) + ( RAYY - ry ) * cos( hudu ) + ry;
	return yr;
}
double BetaRotaryAngleX( double BRAXX, double BRAXY )
{
	double bxr;

	bxr = ( BRAXX - betax ) * cos( betahudu ) - ( BRAXY - betay ) * sin( betahudu ) + betax;
	return bxr;
}
double BetaRotaryAngleY( double BRAYX, double BRAYY )
{
	double byr;
	byr = ( BRAYX - betax ) * sin( betahudu ) + ( BRAYY - betay ) * cos( betahudu ) + betay;
	return byr;
}
void getpot(int NNP,int k,double dblet[200],double source[200]) //Page 10
{
  int j;
  double tx,ty,xn1,yn1,xn2,yn2,xc,yc,pl;
  double l1,l2,twr1,twr2,h,theta;
  xn1=Xnodp[k];
  yn1=Ynodp[k];
  xn2=Xnodp[k+1];
  yn2=Ynodp[k+1];
  pl=Plnth[k];
  tx=(xn2-xn1)/pl;
  ty=(yn2-yn1)/pl;
  for(j=1;j<=NNP;j++){
    xc=Xcolp[j];
    yc=Ycolp[j];
    l1=(xc-xn1)*tx+(yc-yn1)*ty;
    l2=-((xc-xn2)*tx+(yc-yn2)*ty);
    h=-((xc-xn1)*ty-(yc-yn1)*tx);
    twr1=(xc-xn1)*(xc-xn1)+(yc-yn1)*(yc-yn1); /*squre(r1)*/
    twr2=(xc-xn2)*(xc-xn2)+(yc-yn2)*(yc-yn2); /*squre(r2)*/
    theta=atan2(fabs(h)*pl,twr1-pl*l1); //Dblet
    if (j==k){
      dblet[j]=0.5;
      source[j]=0.0;
    }
    else
      dblet[j]=-h/fabs(h)*theta/(2.0*pi);
      source[j]=0.5*((pl-l1)*log(twr2)+l1*log(twr1)-2.0*pl+2.0*fabs(h)*theta)/(2.0*pi);
  }
}

void matrixsolution( double A[ 200 ][ 200 ], double known[ 200 ], int num_var )
{
	int num_eq;
	int k, l;
	int row, column, pivot_column, max_index;
	double max_value, ftemp1, ftemp2, pivot_value;
	double coeff[ 200 ][ 200 ], inv_coeff[ 200 ][ 200 ], var[ 200 ];
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
				if(pivot_column != max_index )
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

