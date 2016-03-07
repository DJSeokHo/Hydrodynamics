#include <stdio.h>
#include <math.h>

#define M 31				//number of surfaces
#define N 200				//Array size
#define PI 3.1415926535		//PI
#define r 0.501				//R of point J
#define a 0.5				//R of cycle
#define GAM 1.0				//GAM

void draw(void);			//mian function

float init_x(float, float);                     //initialization  X of the circle
float init_y(float, float);                     //initialization  Y of the circle
float init_cp(float, float);                    //initialization  Control points of the circle
float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);						//trans radian to degree

float cacu_K_Sm_Sn(float, float, float, float, float, float, float);	//cacu K(Sm, Sn) for Winf
float cacu_rhsm(float, float, float, float);					//cacu rhsm Winf only

void matrix_inversion_using_elementary_operation( float original[N][N], float inverse[N][N], int dim );		//Output GAM(s)

float cacu_rm(float, float, float, float);      //Get Rm
float cacu_Um(float, float, float);             //Get Um
float cacu_Vm(float, float, float);             //Get Vm

int Sorting_min(float R[], int L);

main()
{
	draw();
	return 0;
}

void draw(void)
{
	int i, j, nearest;						//for Loop
	FILE *circle=fopen("circle.dat","wt");			//for geometry
	FILE *pivotalP=fopen("pivotalP.dat","wt");		//for geometry
	FILE *sub=fopen("sub.dat","wt");				//for geometry
						
	float cycx[N]={0,}, cycy[N]={0,}, xmid[N]={0,}, ymid[N]={0,}, deltaSn[N]={0,}, angle, cosbetam[N]={0,}, sinbetam[N]={0,};		//init the cycle

	float Coeff[N]={0,};																//for Martensen
	float K_Sm_Sn[N][N]={0,}, inv_Kmn[N][N]={0,}, Surface_gam[N]={0,}, rhsM[N]={0,};	//for Martensen
	float rhsp[N]={0,}, SIGMA_rhsn[N]={0,};

	float xj, yj;
	float rmj[N]={0,}, Umj[N]={0,}, Vmj[N]={0,};
	
	float temp_rmj;
//*************************************||   Input profile and other data preparation||*************************************
	
	xj = 0.0;
	yj = r;


	angle = 0.0;												//init angle of circle

	for(i = 0; i <= M-1; i++)														//////////////
	{																				//          //
		cycx[i] = init_x(ang_to_rad(angle),a);										//init cycle//
		cycy[i] = init_y(ang_to_rad(angle),a);										//          //
		fprintf(circle, "%0.4f  %0.4f\n", cycx[i], cycy[i]);						//          //
		angle += 360.0 / (M-1);														//////////////
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
		cosbetam[i] = ((cycx[i] - cycx[i+1]))/deltaSn[i];
		sinbetam[i] = ((cycy[i] - cycy[i+1]))/deltaSn[i];
	}	
//*************************************||           Caculate the Umj & Vmj          ||*************************************

	for(i=0;i<M-1;i++)
	{
		rmj[i] = cacu_rm(xmid[i], xj, ymid[i], yj);
		Umj[i] = cacu_Um(ymid[i], yj, rmj[i]);
		Vmj[i] = cacu_Vm(xmid[i], xj, rmj[i]);
	}

//*************************************||        Martensen analysis Winf only       ||*************************************

	for(i=0;i<M-1;i++)																//Martensen analysis
	{
		for(j=0;j<M-1;j++)
		{
			K_Sm_Sn[i][j] = cacu_K_Sm_Sn(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);	//(XM, YM, XN, YN, COS, SIN)
		}
		rhsM[i] = cacu_rhsm(Umj[i], cosbetam[i], Vmj[i], sinbetam[i]);	//(U, cos, V, sin)

	//	printf("%0.4f\n", rhsM[i]);
	}

//*************************************||       rhsm  Back diagonal correction      ||*************************************
	
	nearest = Sorting_min(rmj,M-1);
	rhsM[nearest] = 0.0;
	for(i=0;i<M-1;i++)																//Martensen analysis
	{
		if(i!=nearest)
		{
			rhsM[nearest] += rhsM[i] * deltaSn[i];
		}
	}
	rhsM[nearest] = (-1/deltaSn[nearest]) * rhsM[nearest];
//	printf("%0.4f\n", rhsM[nearest]);

//*************************************||                Matrix inversion           ||*************************************

	matrix_inversion_using_elementary_operation(K_Sm_Sn,inv_Kmn,M-1);
	for(i=0;i<M-1;i++)
	{
		Surface_gam[i]=0.0;
		for(j=0;j<M-1;j++)
		{
			Surface_gam[i]+=inv_Kmn[i][j]*rhsM[j];
		}
		printf("%d -- %0.4f\n", i+1, Surface_gam[i]);
	}

}

float init_x(float cta, float R)
{
	float temp_x,temp_r;
	temp_r=R;
	temp_x=temp_r*(float)cos(cta);
	return temp_x;
}

float init_y(float cta, float R)
{
	float temp_y,temp_r;
	temp_r=R;
	temp_y=temp_r*(float)sin(cta);
	return temp_y;
}

float init_cp(float V_1, float V)
{
	float temp_cp_, temp_v1_, temp_v_;
	temp_v1_ = V_1; temp_v_ = V;
	temp_cp_=0.5*(temp_v1_ + temp_v_);
	return temp_cp_;
}

float ang_to_rad(float ANG)
{
	float temp_rad, temp_ang;
	temp_ang=ANG;
	temp_rad=(temp_ang*PI)/180;
	return temp_rad;
}

float rad_to_ang(float RAD)
{
	float temp_ang_, temp_rad_;
	temp_rad_=RAD;
	temp_ang_ = temp_rad_/PI*180;
	return temp_ang_;
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
		temp_Kmn=COE*(((temp_ym-temp_yn)*temp_cos-(temp_xm-temp_xn)*temp_sin) / ((float)pow((temp_xm-temp_xn),2)+(float)pow((temp_ym-temp_yn),2)));
	}
	
	return temp_Kmn;
}

float cacu_rhsm(float U, float COS, float V, float SIN)
{
	float temp_rhsm, temp_u, temp_cos, temp_v, temp_sin;
	temp_u=U; temp_cos=COS; temp_v=V; temp_sin=SIN;
	temp_rhsm=(-1)*GAM*((temp_u*temp_cos)+(temp_v*temp_sin));
	return temp_rhsm;
}

float cacu_rhsp(float DELTASP, float SIGMA)
{
	float temp_rhsp, temp_deltaSp, temp_sig;
	temp_deltaSp=DELTASP; temp_sig=SIGMA;
	temp_rhsp = (-1)*(1/temp_deltaSp)*temp_sig;
	return temp_rhsp;
}

float cacu_rm(float XM, float X_, float YM, float Y_)
{
	float temp_rm_, temp_xm_, temp_x_, temp_ym_, temp_y_;
	temp_xm_=XM; temp_x_=X_; temp_ym_=YM; temp_y_=Y_;
	temp_rm_ = (float)sqrt( (float)pow((temp_xm_-temp_x_),2) + (float)pow((temp_ym_-temp_y_),2) );
	return temp_rm_;
}

float cacu_Um(float YM, float Y_, float RM_)
{
	float temp_Um_, temp_ym_, temp_y_, temp_rm_;
	temp_ym_=YM; temp_y_=Y_; temp_rm_=RM_;

	if(temp_rm_==0)
	{
		return 0;
	}
	else
	{
		temp_Um_ = ( 1/(2*PI) ) * ( (temp_ym_-temp_y_) / (float)pow(temp_rm_,2) );
		return temp_Um_;
	}
}

float cacu_Vm(float XM, float X_, float RM_)
{
	float temp_Vm_, temp_xm_, temp_x_, temp_rm_;
	temp_xm_=XM; temp_x_=X_; temp_rm_=RM_;

	if(temp_rm_==0)
	{
		return 0;
	}
	else
	{
		temp_Vm_ = ( -1/(2*PI) ) * ( (temp_xm_-temp_x_) / (float)pow(temp_rm_,2) );
		return temp_Vm_;
	}
}

int Sorting_min(float R[], int L)
{
	int i;
	float min=R[0],minN=0;
	for(i=1;i<L;i++)
	{
		if(min>R[i])
		{
			min=R[i];
			minN=i;
		}
	}
	return minN;
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
