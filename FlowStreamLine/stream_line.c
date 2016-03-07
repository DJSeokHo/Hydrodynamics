//*********************************************
//**************   streamline   ***************
//**************      V2.0      ***************
//*********************************************

	

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#define r 1
#define d 2 * r

#define PI 3.141592654
#define M 2 * PI
#define Uinf 1


double cyc( double );
double rzero( double );
double sinceta( double );
double streamline( double, double, double );
void matrix_inversion_using_element_operation( float coeff[20][20], float inv_coeff[20][20], int num_eq );
double CubicEquationssolution( double, double, double );
double realsin( double );
double realcos( double );
double ans( double, double, double );
double QuadraticEquations( double, double, double );

main()
{
	FILE *drawstreaml = fopen( "streaml.dat", "wt" );

	int i, j;
	double CEX, CEYUP, CEYDOWN;
	double rz, psi;

	double cyup, cydown;
	double x;
	/***************************************************///drawcylinder
	printf( "****************************************************\n" );
	printf( "**************  draw cylinder  *********************\n" );
	printf( "****************************************************\n" );
	
	fprintf( drawstreaml, "\nTitle=“cylinder”\n");
	fprintf( drawstreaml, "Variables=“X”,“Y”,“Z”\n");
	fprintf( drawstreaml, "Zone I = %d \n\n", 20 * d + 1 );

	for( i = -20 * r; i <= 20 * r; i++ ){
		x = ( double ) i / 20;
		cyup = cyc( x );
		cydown = cyup * -1;

//		printf( "%0.5lf", x );
//		printf( "   %0.5lf   %0.5lf\n", cyup, cydown );

		fprintf( drawstreaml, "%0.5lf", x );
		fprintf( drawstreaml, "   %0.5lf    %0.5lf\n", cyup, cydown );
	}

	/************************************************///streamfuc
	printf( "****************************************************\n" );
	printf( "****************  draw stream lines  ***************\n" );
	printf( "****************************************************\n" );

	for( i = 1; i <= 8; i++ )
	{
		rz = rzero( 2 );  //n=2
	
		psi = ( double ) i / 4;

		fprintf( drawstreaml, "\nTitle=“streaml”\n");
		fprintf( drawstreaml, "Variables=“X”,“Y”,“Z”\n");
		fprintf( drawstreaml, "Zone I=120\n\n");
		for( j = -60; j <= 60; j++ ){
			CEX = ( double ) j / 10;
			if( CEX == 0 )
			{
				CEYUP = QuadraticEquations( Uinf, psi, rz );
				CEYDOWN = CEYUP * -1;
//				printf( "%0.5lf", CEX );
//				printf( "   %0.5lf   %0.5lf\n",CEYUP, CEYDOWN );
				fprintf( drawstreaml, "%0.5lf   %0.5lf   %0.5lf\n",CEX, CEYUP, CEYDOWN );
			}
			else{
				CEYUP = CubicEquationssolution( CEX, rz, psi );
				CEYDOWN = CEYUP * -1;
//				printf( "%0.5lf", CEX );
//				printf( "   %0.5lf   %0.5lf\n",CEYUP, CEYDOWN );
				fprintf( drawstreaml, "%0.5lf   %0.5lf   %0.5lf\n",CEX, CEYUP, CEYDOWN );
			}
		}
	}
	fclose( drawstreaml );
		/************************************************///potentialfuc

	return 0;
}

double rzero( double n )
{
	double temp = 0;
	double rzero = 0;

	temp = n * PI * Uinf;
	rzero = M / temp;
	rzero = sqrt( rzero );
	
	return rzero;
}

double sinceta( double n )
{
	double temp = 0;
	double sc = 0;

	return 0;
}

double cyc( double x )
{
	double cr = 0, cx = 0, temp = 0;
	cr = pow( r, 2 );
	cx = pow( x, 2 );
	temp = cr - cx;
	temp = sqrt( temp );
	return temp;
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

double CubicEquationssolution( double x, double rz, double psi )
{
	double xz = 0, fz = 0, fo = 0, ft = 0;
	double xo = -100, xt = 100;

	do{
		fo = ( Uinf * ( -1 ) ) * xo * xo * xo + psi * xo * xo + ( ( Uinf * x * x - rz ) * ( -1 ) ) * xo + ( psi * x * x ); //x,rz: variable
		if( fo == 0 )
		{
			xz = xo;
			goto loop;
		}
		ft = ( Uinf * ( -1 ) ) * xt * xt * xt + psi * xt * xt + ( ( Uinf * x * x - rz ) * ( -1 ) ) * xt + ( psi * x * x ); //x,rz: variable
		if( ft == 0 )
		{
			xz = xt;
			goto loop;
		}
	}while( fo * ft > 0 );
	do{
		xz = ( xo + xt ) / 2;
		fz = ( Uinf * ( -1 ) ) * xz * xz * xz + psi * xz * xz + ( ( Uinf * x * x - rz ) * ( -1 ) ) * xz + ( psi * x * x );
		if( fz == 0 )
			break;
		if( fz * fo < 0 )
		{
			xt = xz;
			ft = fz;
		}
		else
		{
			xo = xz;
			fo = fz;
		}
	}while( fabs( fz ) >= 1e-5 );

	
loop:	return xz;
}

double ans(double ansa,double ansb,double ansc)
{
	double ansd, ansx1, ansx2, anss;
	ansd = ansb * ansb - 4 * ansa * ansc;
	if( ansd < 0 )
	{
		printf( "no solution!\n" );
		return 0;
	}
	else
	{
		anss = sqrt( ansd );
		ansx1 = ( -ansb + anss ) / ( 2 * ansa ), ansx2 = ( -ansb -anss ) / ( 2 * ansa );
		return ansx1;
	}
}

double QuadraticEquations( double QEa, double QEb, double QEc )
{
	double QEX;
	QEa = QEa * 1;
	QEb = QEb * -1;
	QEc = QEc * -1;
	if( QEa != 0 )
	{
		QEX = ans( QEa, QEb, QEc );
		return QEX;
	}
	else
	{
		printf( "no solution!" );
		return 0;
	}
}

void matrix_inversion_using_element_operation( float coeff[20][20], float inv_coeff[20][20], int num_eq )
{
	{
		int row, column, pivot_column, max_index;
		float max_value, ftemp1, ftemp2, pivot_value;

		for( row = 0; row < num_eq; row++ )
			for( column = 0; column < num_eq; column++ ){
				if( row == column )
					inv_coeff[row][column] = 1;
				else
					inv_coeff[row][column] = 0;
			}

			for( pivot_column = 0; pivot_column < num_eq; pivot_column++ ){
				max_index = coeff[0][column];
				max_value = 0;

				for( row = pivot_column; row < num_eq; row++ )
					if( coeff[row][pivot_column]*coeff[row][pivot_column] > max_value*max_value ){
						max_index = row;
						max_value = coeff[row][pivot_column];
					}

					if(pivot_column != max_index )
						for( column = 0; column < num_eq; column++ ){
							ftemp1 = coeff[pivot_column][column];
							ftemp2 = inv_coeff[pivot_column][column];

							coeff[pivot_column][column] = coeff[max_index][column];
							
							inv_coeff[pivot_column][column] = inv_coeff[max_index][column];
							coeff[max_index][column] = ftemp1;
							inv_coeff[max_index][column] = ftemp2;
						}

						pivot_value = coeff[pivot_column][pivot_column];
						for(column = 0; column < num_eq; column++ ){
							coeff[pivot_column][column] /= pivot_value;
							inv_coeff[pivot_column][column] /= pivot_value;

						}

						for( row = 0; row < num_eq; row++ )
							if( row != pivot_column ){
								ftemp1 = coeff[row][pivot_column];
								for( column = 0; column < num_eq; column++ ){

									coeff[row][column] -= ftemp1*coeff[pivot_column][column];
									inv_coeff[row][column] -= ftemp1*inv_coeff[pivot_column][column];

								}
							}
			}
	}
}
