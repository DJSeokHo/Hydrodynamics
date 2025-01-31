#include <stdio.h>
#include <math.h>

#define N 19
#define PI 3.1415926535
#define a 0.5

void matrix_inversion_using_elementary_operation( float original[N][N], float inverse[N][N], int dim );
float cacu_beta(float,float,float,float,float,float);	//Get angle beta

float init_x(float, float);                     //initialization  X of the circle
float init_y(float, float);                     //initialization  Y of the circle
float init_cp(float, float);                    //initialization  Control points of the circle
float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);

float cacu_K_Sm_Sn(float, float, float, float, float, float, float);	//cacu K(Sm, Sn), (XM, YM, XN, YN, COS, SIN)
float cacu_rhsm(float, float, float, float);					//cacu rhsm
float cacu_slope(float, float, float, float);		//cacu betam(x1,x2,y1,y2)

main()
{
	int i, j;

	FILE *circle=fopen("circle.dat","wt");
	FILE *pivotalP=fopen("pivotalP.dat","wt");

	float Winf, Uinf, Vinf, flow_angle;															//init flow

	float cycx[N]={0,}, cycy[N]={0,}, xmid[N]={0,}, ymid[N]={0,}, deltaSm[N]={0,}, angle;				//init the cycle

	float K_Sm_Sn[N][N]={0,}, inv_Kmn[N][N]={0,}, Surface_gam[N]={0,}, rhsM[N]={0,}, slope[N]={0,}, cosbetam[N]={0,}, sinbetam[N]={0,};//init surface gam


//--------------------------------------------------------------------------------------------------------
//**********************************||       init Uniform flow       ||***********************************
	flow_angle=0;
	Winf = 1.0;
	Uinf = Winf*(float)cos(ang_to_rad(flow_angle));
	Vinf = Winf*(float)sin(ang_to_rad(flow_angle));
//********************************************************************************************************


//--------------------------------------------------------------------------------------------------------
//***********************************||       init the cycle       ||************************************
	for(i = 0, angle = 0; i <= N-1, angle <= 360; i++, angle += (360/(N-1)))
	{
		cycx[i] = init_x(ang_to_rad(angle),a);
		cycy[i] = init_y(ang_to_rad(angle),a);
		fprintf(circle, "%0.4f  %0.4f\n", cycx[i], cycy[i]);
	}
//********************************************************************************************************


//--------------------------------------------------------------------------------------------------------
//****************************||         init cosbeta & sinbeta         ||********************************
	for(i = 0; i < N - 1; i++)
	{
		slope[i] = cacu_slope(cycx[i+1], cycx[i], cycy[i+1], cycy[i]);
//		printf("%d -- %0.4f\n", i+1, slope[i]);
		cosbetam[i] = (float)cos(slope[i]);
		sinbetam[i] = (float)sin(slope[i]);
	}
//********************************************************************************************************


//--------------------------------------------------------------------------------------------------------
//*****************************||         init mid point and Vs        ||*********************************
	for(i = 0; i < N - 1; i++)
	{
		xmid[i] = init_cp(cycx[i+1],cycx[i]);
		ymid[i] = init_cp(cycy[i+1],cycy[i]);
		fprintf(pivotalP, "%0.4f  %0.4f\n", xmid[i], ymid[i]);
	}
//********************************************************************************************************


//--------------------------------------------------------------------------------------------------------
//*************************||           cacu element length           ||**********************************
	for(i = 0; i < N - 1; i++)	//element's deltaS
	{
		deltaSm[i] = (float)sqrt((float)pow((cycx[i+1] - cycx[i]),2) + (float)pow((cycy[i+1] - cycy[i]),2));
	}
//********************************************************************************************************

//--------------------------------------------------------------------------------------------------------
//**************************||           init the K(Sm,Sn)           ||***********************************
	for(i=0;i<N-1;i++)
	{
		for(j=0;j<N-1;j++)
		{
			K_Sm_Sn[i][j] = cacu_K_Sm_Sn(xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);// * deltaSm[j];	//(XM, YM, XN, YN, COS, SIN)
		}
	}
//********************************************************************************************************


//--------------------------------------------------------------------------------------------------------
//****************************||           init the rhsm           ||*************************************
	for(i=0;i<N-1;i++)
	{
		rhsM[i] = cacu_rhsm(Uinf, cosbetam[i], Vinf, sinbetam[i]);	//(U, cos, V, sin)
	}
//********************************************************************************************************


//--------------------------------------------------------------------------------------------------------
//*************************||           cacu the surface GAM           ||*********************************
	matrix_inversion_using_elementary_operation(K_Sm_Sn,inv_Kmn,N-1);

	for(i=0;i<N-1;i++)
	{
		Surface_gam[i]=0.0;
		for(j=0;j<N-1;j++)
			Surface_gam[i]+=inv_Kmn[i][j]*rhsM[j];
	//	printf("%0.4f\n", Surface_gam[i]);
	}
//********************************************************************************************************



	return 0;
}

float init_x(float cta, float R)
{
	float temp_x,temp_r;
	temp_r=R;
	temp_x=temp_r*(1-(float)cos(cta));
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

float cacu_K_Sm_Sn(float COE, float XM, float YM, float XN, float YN, float COS, float SIN)
{
	float temp_Kmn, temp_xm, temp_ym, temp_xn, temp_yn, temp_cos, temp_sin, temp_coe;
	temp_xm=XM; temp_ym=YM; temp_xn=XN; temp_yn=YN; temp_cos=COS; temp_sin=SIN; temp_coe=COE;
	if((temp_xm==temp_xn) && (temp_ym==temp_yn))
	{
		temp_Kmn=-0.5;
	}
	else
	{
		temp_Kmn=COE*(((temp_ym-temp_yn)*temp_cos-(temp_xm-temp_xn)*temp_sin)/((float)pow((temp_xm-temp_xn),2)+(float)pow((temp_ym-temp_yn),2)));
	}
	return temp_Kmn;
	temp_Kmn=0; temp_xm=0; temp_ym=0; temp_xn=0; temp_yn=0; temp_cos=0; temp_sin=0; temp_coe=0;
}

float cacu_rhsm(float U, float COS, float V, float SIN)
{
	float temp_rhs, temp_u, temp_cos, temp_v, temp_sin;
	temp_u=U; temp_cos=COS; temp_v=V; temp_sin=SIN;
	temp_rhs=(temp_u*temp_cos)-(temp_v*temp_sin);
	return temp_rhs;
	temp_rhs=0; temp_u=0; temp_cos=0; temp_v=0; temp_sin=0;
}

float cacu_slope(float X1, float X2, float Y1, float Y2)
{
	float temp_slope, temp_x1, temp_x2, temp_y1, temp_y2;
	temp_x1=X1; temp_x2=X2; temp_y1=Y1; temp_y2=Y2;
	temp_slope=(temp_y2-temp_y1)/(temp_x2-temp_x1);
	return temp_slope;
	temp_slope=0; temp_x1=0; temp_x2=0; temp_y1=0; temp_y2=0;
}

float cacu_beta(float X1, float Y1, float X0, float Y0, float X2, float Y2)
{
	float temp_x1_, temp_y1_, temp_x0_, temp_y0_, temp_x2_, temp_y2_;
	float temp_rad_, temp_beta_;
	temp_x1_=X1; temp_y1_=Y1; temp_x0_=X0; temp_y0_=Y0; temp_x2_=X2; temp_y2_=Y2;
	temp_rad_ = (float)fabs((float)atan((temp_y0_-temp_y1_)/(temp_x0_-temp_x1_))-(float)atan((temp_y2_-temp_y0_)/(temp_x2_-temp_x0_)));
	temp_beta_ = temp_rad_ * (180 / PI);
	if(((temp_x1_<0) && (temp_x2_>0))||((temp_x1_>0) && (temp_x2_<0)))
		return 180 - temp_beta_;
	else
		return temp_beta_;
	temp_x1_=0; temp_y1_=0; temp_x0_=0; temp_y0_=0; temp_x2_=0; temp_y2_=0;
	temp_rad_=0, temp_beta_=0;
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
