#include <stdio.h>
#include <math.h>

#define M 11				//number of surfaces
#define N 200				//Array size
#define PI 3.1415926535		//PI
#define Winf 1.0			//Winf
#define a 0.5				//R of cycle
#define GAM 1.0				//GAM
#define flow_angle 0.0		//flow angle

void draw(void);			//mian function

float init_x(float, float);                     //initialization  X of the circle
float init_y(float, float);                     //initialization  Y of the circle
float init_cp(float, float);                    //initialization  Control points of the circle
float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);						//trans radian to degree
float cacu_slope(float, float, float, float);		//cacu betam(x1,x2,y1,y2)

float cacu_K_Sm_Sn_Winf(float, float, float, float, float, float, float);	//cacu K(Sm, Sn) for Winf
float cacu_rhsm_Winf(float, float, float, float);					//cacu rhsm Winf only

float cacu_K_Sm_Sn_ALL(float, float, float, float, float, float, float);	//cacu K(Sm, Sn) for All
float cacu_rhsm_ALL(float, float, float, float);					//cacu rhsm	all 

void Crout_method(float original[N][N], float L[N][N], float U[N][N], int dim);						//Output GAM(s)
void matrix_inversion_using_elementary_operation( float original[N][N], float inverse[N][N], int dim );		//Output GAM(s)

main()
{
	draw();
	return 0;
}

void draw(void)
{
	int i, j, k, time, T;								//for Loop
	FILE *circle=fopen("circle.dat","wt");			//for geometry
	FILE *pivotalP=fopen("pivotalP.dat","wt");		//for geometry
	FILE *shedP=fopen("shedP.dat","wt");			//for geometry

	float Uinf, Vinf;								//init Winf
	float cycx[N]={0,}, cycy[N]={0,}, xmid[N]={0,}, ymid[N]={0,}, deltaSn[N]={0,}, angle, slope[N]={0,}, cosbetam[N]={0,}, sinbetam[N]={0,};		//init the cycle
	float epsl, deltaGAM_A_angle, deltaGAM_A_x, deltaGAM_A_y, deltaGAM_A;				//init deltaGAM A

	float Coeff[N]={0,};																//for Martensen
	float K_Sm_Sn_Winf[N][N]={0,}, inv_Kmn_Winf[N][N]={0,}, Surface_gam_Winf[N]={0,}, rhsM_Winf[N]={0,};	//for Martensen
	float L[N][N], U[N][N], ftemp, ytemp[N];

	float deltaT=0.15;																	//delta T

//*************************************||   Input profile and other data preparation||*************************************

	Uinf = Winf*(float)cos(ang_to_rad(flow_angle));				//init Uniform flow
	Vinf = Winf*(float)sin(ang_to_rad(flow_angle));				//init Uniform flow

	for(i = 0, angle = 0; i <= M-1, angle <= 360; i++, angle += (360/(M-1)))		//////////////
	{																				//          //
		cycx[i] = init_x(ang_to_rad(angle),a);										//init cycle//
		cycy[i] = init_y(ang_to_rad(angle),a);										//          //
		fprintf(circle, "%0.4f  %0.4f\n", cycx[i], cycy[i]);						//          //
	}																				//////////////

	for(i = 0; i < M - 1; i++)														//////////////
	{																				//          //
		xmid[i] = init_cp(cycx[i+1],cycx[i]);										// init mid //
		ymid[i] = init_cp(cycy[i+1],cycy[i]);										//   point  //
		fprintf(pivotalP, "%0.4f  %0.4f\n", xmid[i], ymid[i]);						//          //
	}																				//////////////																			//////////////

	for(i = 0; i < M - 1; i++)														//element's deltaS
	{																				//Coefficient of Martensen
		deltaSn[i] = (float)sqrt((float)pow((cycx[i+1] - cycx[i]),2) + (float)pow((cycy[i+1] - cycy[i]),2));
		Coeff[i] = deltaSn[i]/(2*PI);
	}
	
	for(i = 0; i < M - 1; i++)														//////////////
	{																				//   init   //										// sinbetam //
		cosbetam[i] = ymid[i]/deltaSn[i];											//////////////
		sinbetam[i] = xmid[i]/deltaSn[i];			// In here function cos&sin can only take a rad but not a slope!!!
	}	

//*************************************||        Martensen analysis Winf only       ||*************************************

	for(i=0;i<M-1;i++)																//Martensen analysis
	{
		for(j=0;j<M-1;j++)
		{
			K_Sm_Sn_Winf[i][j] = cacu_K_Sm_Sn_Winf(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);	//(XM, YM, XN, YN, COS, SIN)
		}
		rhsM_Winf[i] = cacu_rhsm_Winf(Uinf, cosbetam[i], Vinf, sinbetam[i]);	//(U, cos, V, sin)
	}

	matrix_inversion_using_elementary_operation(K_Sm_Sn_Winf,inv_Kmn_Winf,M-1);
	for(i=0;i<M-1;i++)
	{
		Surface_gam_Winf[i]=0.0;
		for(j=0;j<M-1;j++)
			Surface_gam_Winf[i]+=inv_Kmn_Winf[i][j]*rhsM_Winf[j];
		printf("%d -- %0.4f\n", i+1, Surface_gam_Winf[i]);
	}

//*************************************************************************************************************************
	for(time=0;time<1;time++)
	{
		T = 1;
		epsl = deltaSn[0]/2;
		deltaGAM_A_angle = ang_to_rad(45);
		deltaGAM_A_x = epsl * (float)cos(deltaGAM_A_angle) + xmid[2];
		deltaGAM_A_y = epsl * (float)sin(deltaGAM_A_angle) + ymid[2];
		fprintf(shedP, "%0.4f %0.4f\n", deltaGAM_A_x, deltaGAM_A_y);
		deltaGAM_A = 0.5 * (float)pow(Surface_gam_Winf[3],2) * deltaT;
		printf("%0.4f\n", deltaGAM_A);
		//***********************||   Shed two free vortex elements deltaGAM A, deltaGAM B   ||****************************
		for(i=0;i<(M-1)+T;i++)
		{
			//*************************||   Martensen analysis Winf and all shed vortices   ||*****************************
			
			//**************************************||     Vortex convection     ||****************************************
		}

		//*****************************************||   Advance delta T   ||***********************************************
	}
}




float init_x(float cta, float R)
{
	float temp_x,temp_r;
	temp_r=R;
//	temp_x=temp_r*(1-(float)cos(cta));
	temp_x=temp_r*(float)cos(cta);
	return temp_x;
	temp_r=0; temp_x=0;
}

float init_y(float cta, float R)
{
	float temp_y,temp_r;
	temp_r=R;
	temp_y=temp_r*(float)sin(cta);
	return temp_y;
	temp_r=0; temp_y=0;
}

float init_cp(float V_1, float V)
{
	float temp_cp_, temp_v1_, temp_v_;
	temp_v1_ = V_1; temp_v_ = V;
	temp_cp_=0.5*(temp_v1_ + temp_v_);
	return temp_cp_;
	temp_cp_ = 0; temp_v1_ = 0; temp_v_ = 0;
}

float ang_to_rad(float ANG)
{
	float temp_rad, temp_ang;
	temp_ang=ANG;
	temp_rad=(temp_ang*PI)/180;
	return temp_rad;
	temp_rad=0; temp_ang=0;
}

float rad_to_ang(float RAD)
{
	float temp_ang_, temp_rad_;
	temp_rad_=RAD;
	temp_ang_ = temp_rad_/PI*180;
	return temp_ang_;
	temp_rad_=0; temp_ang_=0;
}

float cacu_K_Sm_Sn_Winf(float COE, float XM, float YM, float XN, float YN, float COS, float SIN)
{
	float temp_Kmn, temp_xm, temp_ym, temp_xn, temp_yn, temp_cos, temp_sin, temp_coe;
	temp_xm=XM; temp_ym=YM; temp_xn=XN; temp_yn=YN; temp_cos=COS; temp_sin=SIN; temp_coe=COE;
	if((temp_xm==temp_xn) && (temp_ym==temp_yn))
	{
		temp_Kmn=-0.5;
	}
	else
	{
		temp_Kmn=COE*(((temp_ym-temp_yn)*temp_cos-(temp_xm-temp_xn)*temp_sin) / ((float)pow((temp_xm-temp_xn),2)+(float)pow((temp_ym-temp_yn),2)));
	}
	
	return temp_Kmn;
	temp_Kmn=0; temp_xm=0; temp_ym=0; temp_xn=0; temp_yn=0; temp_cos=0; temp_sin=0; temp_coe=0;
}

float cacu_rhsm_Winf(float U, float COS, float V, float SIN)
{
	float temp_rhs, temp_u, temp_cos, temp_v, temp_sin;
	temp_u=U; temp_cos=COS; temp_v=V; temp_sin=SIN;
	temp_rhs=((-1)*(temp_u*temp_cos))+((-1)*(temp_v*temp_sin));
	return temp_rhs;
	temp_rhs=0; temp_u=0; temp_cos=0; temp_v=0; temp_sin=0;
}

void matrix_inversion_using_elementary_operation( float original[N][N], float inverse[N][N], int dim )
{
	int row, column, pivot_column, max_index;
	float max_value, ftemp1, ftemp2, pivot_value;
	for( row = 0; row < dim; row++ )
		for( column = 0; column < dim; column++ ){
			if( row == column )
				inverse[row][column] = 1;
			else
				inverse[row][column] = 0;
		}
		for( pivot_column = 0; pivot_column < dim; pivot_column++ ){
			max_index = original[0][column];
			max_value = 0;
			for( row = pivot_column; row < dim; row++ )
				if( original[row][pivot_column]*original[row][pivot_column] > max_value*max_value ){
					max_index = row;
					max_value = original[row][pivot_column];
				}
				if(pivot_column != max_index )
					for( column = 0; column < dim; column++ ){
						ftemp1 = original[pivot_column][column];
						ftemp2 = inverse[pivot_column][column];
						original[pivot_column][column] = original[max_index][column];
						inverse[pivot_column][column] = inverse[max_index][column];
						original[max_index][column] = ftemp1;
						inverse[max_index][column] = ftemp2;
					}
					pivot_value = original[pivot_column][pivot_column];
					for(column = 0; column < dim; column++ ){
						original[pivot_column][column] /= pivot_value;
						inverse[pivot_column][column] /= pivot_value;
					}
					for( row = 0; row < dim; row++ )
						if( row != pivot_column ){
							ftemp1 = original[row][pivot_column];
							for( column = 0; column < dim; column++ ){
								original[row][column] -= ftemp1*original[pivot_column][column];
								inverse[row][column] -= ftemp1*inverse[pivot_column][column];
							}
						}
		}
}

void Crout_method(float original[N][N], float L[N][N], float U[N][N], int dim)
{
		int i, row, column;
	float ftemp;
	for(row=0;row<dim;row++)
		L[row][0]=original[row][0];
	for(row=0;row<dim;row++)
		U[row][row]=1.0;
	for(column=1;column<dim;column++)
		U[0][column]=original[0][column]/original[0][0];
	for(column=1;column<dim;column++)
	{
		for(row=column;row<dim;row++)
		{
			ftemp=0.0;
			for(i=0;i<column;i++)
				ftemp+=L[row][i]*U[i][column];
			L[row][column]=original[row][column]-ftemp;
		}
		for(row=column+1;row<dim;row++)
		{
			ftemp=0;
			for(i=0;i<column;i++)
				ftemp+=L[column][i]*U[i][row];
			ftemp=original[column][row]-ftemp;
			U[column][row]=ftemp/L[column][column];
		}
	}
}
