#include <stdio.h>
#include <math.h>

#define M 9
#define N 10
#define PI 3.1415926535
#define a 0.5
#define Winf 1.0
#define flow_angle 0.0
#define shedsurface 1

float init_x(float, float);                     //initialization  X of the circle
float init_y(float, float);                     //initialization  Y of the circle
float init_cp(float, float);                    //initialization  Control points of the circle
float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);

float cacu_rhsm_Winf(float, float, float, float);					//cacu rhsm
float cacu_K_Sm_Sn(float, float, float, float, float, float, float);	//cacu K(Sm, Sn) for MJ
float cacu_rhsm(float, float, float, float, float);					//cacu rhsm for MJ 

void matrix_inversion_using_elementary_operation_Winf( float original[N][N], float inverse[N][N], int dim );		//Output GAM(s)
void matrix_inversion_using_elementary_operation( float original[N][N], float inverse[N][N], int dim );		//Output GAM(s)

float cacu_rm(float, float, float, float);      //Get Rmj
float cacu_Um(float, float, float);             //Get Umj
float cacu_Vm(float, float, float);             //Get Vmj
int Sorting_min(float R[], int L);

float cacu_Umn(float, float, float, float);		//Get Umn. Here m & n means convection shed points only ?
float cacu_Vmn(float, float, float, float);		//Get Vmn. Here m & n means convection shed points only ?

main()
{
	int i, j, vortexNum, time, nearest[N];

	FILE *circle=fopen("circle.dat","wt");				//for geometry
	FILE *pivotalP=fopen("pivotalP.dat","wt");			//for geometry
	FILE *driftPath=fopen("driftPath.dat","wt");		
	FILE *shed=fopen("shed.dat","wt");

	FILE *mix_A=fopen("martixA.dat","wt");
	FILE *mix_B=fopen("martixB.dat","wt");
	FILE *mix_X=fopen("martixX.dat","wt");

	FILE *velocity=fopen("velocity.dat","wt");
	FILE *vortexGAM=fopen("vortexGAM.dat","wt");

	float Uinf, Vinf;
	float cycx[N]={0,}, cycy[N]={0,}, xmid[N]={0,}, ymid[N]={0,}, deltaSn[N]={0,}, angle;				//init the cycle
	float x[N], y[N], tx[N], ty[N], nx[N], ny[N], cosbetam[N]={0,}, sinbetam[N]={0,};

	float Surface_gam_Winf[N]={0,}, rhsM_Winf[N]={0,};
	float K_Sm_Sn[N][N]={0,}, inv_Kmn[N][N]={0,}, Surface_gam[N]={0,}, rhsM[N]={0,};
	
	float rmj[N]={0,}, Umj[N]={0,}, Vmj[N]={0,};	//for induced-velocity
	float Umn[N]={0,}, Vmn[N]={0,};					//for self-induced-velocity

	float delta_GAM[N], deltaT, offset;				//all shed points strength and delta T

	float shed_x[N]={0,}, shed_y[N]={0,};			//all shed vortices									

	float point_x[N]={0,}, point_y[N]={0,};			//all shed vortices								

	float point_x_next[N]={0,}, point_y_next[N]={0,};							//for drift path

	float sigma_rhsm_for_mj[N]={0,}, sigma_rhsm_mj[N]={0,};					//for K function's second sigma item

	float udm[N]={0,}, vdm[N]={0,};
	float udm_sig_one[N]={0,}, udm_sig_two[N]={0,}, vdm_sig_one[N]={0,}, vdm_sig_two[N]={0,};

	float outputstreamL_x[N]={0,}, outputstreamL_y[N]={0,};

//--------------------------------------------------------------------------------------------------------
//**********************************||       init Uniform flow       ||***********************************
	Uinf = Winf*(float)cos(ang_to_rad(flow_angle));
	Vinf = Winf*(float)sin(ang_to_rad(flow_angle));
	printf("Uinf = %0.4f, Vinf = %0.4f\n", Uinf, Vinf);
//********************************************************************************************************

//--------------------------------------------------------------------------------------------------------
//***********************************||       init the cycle       ||************************************
	angle = 0.0;
	for(i = 0; i <= M-1; i++)														//////////////
	{																				//          //
		cycx[i] = init_x(ang_to_rad(angle),a);										//init cycle//
		cycy[i] = init_y(ang_to_rad(angle),a);										//          //
		fprintf(circle, "%0.4f  %0.4f\n", cycx[i], cycy[i]);						//          //
		angle += 360.0 / (M-1);														//////////////
	}
//********************************************************************************************************

//--------------------------------------------------------------------------------------------------------
//*****************************||         init mid point and Vs        ||*********************************
	for(i = 0; i < M - 1; i++)
	{
		xmid[i] = init_cp(cycx[i+1],cycx[i]);
		ymid[i] = init_cp(cycy[i+1],cycy[i]);
		fprintf(pivotalP, "%0.4f  %0.4f\n", xmid[i], ymid[i]);
	}
//********************************************************************************************************

//--------------------------------------------------------------------------------------------------------
//****************************||         init cosbeta & sinbeta         ||********************************
	for(i = 0; i < M - 1; i++)	//element's deltaS
	{
		x[i] = cycx[i+1] - cycx[i];
		y[i] = cycy[i+1] - cycy[i];
		deltaSn[i] = (float)sqrt((float)pow((cycx[i+1] - cycx[i]),2) + (float)pow((cycy[i+1] - cycy[i]),2));
		tx[i] = x[i] / deltaSn[i];
		ty[i] = y[i] / deltaSn[i];
		nx[i] = ty[i];
		ny[i] = -tx[i];
	}
	for(i = 0; i < M - 1; i++)
	{
		cosbetam[i] = -tx[i];						
		sinbetam[i] = -ty[i];
	//	printf("%d  cosbeta = %0.4f , sinbeta =  %0.4f\n", i+1, cosbetam[i], sinbetam[i]);
	}
//********************************************************************************************************

//--------------------------------------------------------------------------------------------------------
//**************************||           init the K(Sm,Sn)           ||***********************************
	for(i=0;i<M-1;i++)
	{
		for(j=0;j<M-1;j++)
		{
			K_Sm_Sn[i][j] = cacu_K_Sm_Sn(deltaSn[i], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]) + deltaSn[j];	//(XM, YM, XN, YN, COS, SIN)
			fprintf(mix_A, "%0.6f ", K_Sm_Sn[i][j]);
		}
		fprintf(mix_A, "\n");
	}
//********************************************************************************************************

//--------------------------------------------------------------------------------------------------------
//****************************||           init the rhsm           ||*************************************
	for(i=0;i<M-1;i++)
	{
		rhsM_Winf[i] = cacu_rhsm_Winf(Uinf, cosbetam[i], Vinf, sinbetam[i]);	//(U, cos, V, sin)
		fprintf(mix_B, "%0.6f\n", rhsM_Winf[i]);
	}
//********************************************************************************************************


//--------------------------------------------------------------------------------------------------------
//*************************||           cacu the surface GAM           ||*********************************
	matrix_inversion_using_elementary_operation_Winf(K_Sm_Sn,inv_Kmn,M-1);

	for(i=0;i<M-1;i++)
	{
		Surface_gam_Winf[i]=0.0;
		for(j=0;j<M-1;j++)
		{
			Surface_gam_Winf[i]+=inv_Kmn[i][j]*rhsM_Winf[j];
		}
		fprintf(mix_X, "%0.4f\n", Surface_gam_Winf[i]);
	}
//********************************************************************************************************

//*************************************||            Shed all free vortex           ||*************************************
	deltaT = 0.098;
	offset = 0.5;

	time = 1;

	vortexNum = time * shedsurface;

/*	for(i=0;i<M-1;i++)
	{
		shed_x[i] = xmid[i] + nx[i] * offset * deltaSn[i];
		shed_y[i] = ymid[i] + ny[i] * offset * deltaSn[i];
		delta_GAM[i] = Surface_gam_Winf[i] * deltaSn[i];	//Get surface gam's strength
		fprintf(shed, "%0.6f %0.6f\n", shed_x[i], shed_y[i]);
		point_x[i] = shed_x[i];
		point_y[i] = shed_y[i];
	}*/

	shed_x[0] = xmid[1] + nx[1] * offset * deltaSn[1];
	shed_y[0] = ymid[1] + ny[1] * offset * deltaSn[1];
	fprintf(shed, "%0.6f %0.6f\n", shed_x[0], shed_y[0]);
	delta_GAM[0] = Surface_gam_Winf[0] * deltaSn[0];	//Get surface gam's strength
	point_x[0] = shed_x[0];
	point_y[0] = shed_y[0];

	fprintf(mix_A, "\nWinf and all shedding vortices\n\n");
	fprintf(mix_B, "\nWinf and all shedding vortices\n\n");
	fprintf(mix_X, "\nWinf and all shedding vortices\n\n");
//*************************************||           Caculate the Umj & Vmj          ||*************************************
	for(i=0;i<M-1;i++)
	{
		sigma_rhsm_mj[i] = 0.0;
		for(j=0;j<vortexNum;j++)
		{
			rmj[j] = cacu_rm(xmid[i], point_x[j], ymid[i], point_y[j]);
			Umj[j] = cacu_Um(ymid[i], point_y[j], rmj[j]);
			Vmj[j] = cacu_Vm(xmid[i], point_x[j], rmj[j]);
			sigma_rhsm_for_mj[j] = delta_GAM[j] * (1.0 + Umj[j] * cosbetam[i] + Vmj[j] * sinbetam[i]);
			sigma_rhsm_mj[i] += sigma_rhsm_for_mj[j];
			nearest[i] = Sorting_min(rmj, vortexNum);
		}
	}
//********************************************************************************************************


//*************************||           cacu the new surface GAM           ||*********************************
	for(i=0;i<M-1;i++)
	{
		rhsM[i] = rhsM_Winf[i] - sigma_rhsm_mj[i];
		fprintf(mix_B, "%0.6f\n", rhsM[i]);
	}
	fprintf(mix_B, "\n");

	for(i=0;i<M-1;i++)
	{
		Surface_gam[i]=0.0;
		for(j=0;j<M-1;j++)
		{
			Surface_gam[i]+=inv_Kmn[i][j]*rhsM[j];				//Get surface GAM
		}
	//	printf("%0.4f\n", Surface_gam[i]);
		fprintf(mix_X, "%0.6f\n", Surface_gam[i]);
	}
	fprintf(mix_X, "\n");
//********************************************************************************************************

//*************************************||     Cacu the convection of Umn & Vmn      ||*************************************
	for(i=0;i<vortexNum;i++)
	{
		udm_sig_one[i] = 0.0; 
		vdm_sig_one[i] = 0.0;

		for(j=0;j<vortexNum;j++)
		{
			if(i!=j)
			{
				udm_sig_one[i] += delta_GAM[i] * cacu_Umn(point_x[i], point_y[i], point_x[j], point_y[j]);
				vdm_sig_one[i] += delta_GAM[i] * cacu_Vmn(point_x[i], point_y[i], point_x[j], point_y[j]);
			}
		}
	//	printf("%0.4f, %0.4f\n", udm_sig_one[i], vdm_sig_one[i]);
	}
	
	for(i=0;i<vortexNum;i++)
	{
		udm_sig_two[i] = 0.0;
		vdm_sig_two[i] = 0.0;
	
		for(j=0;j<M-1;j++)									//init two item
		{
			rmj[j] = cacu_rm(xmid[i], point_x[j], ymid[i], point_y[j]);
			Umj[j] = cacu_Um(ymid[i], point_y[j], rmj[j]);
			Vmj[j] = cacu_Vm(xmid[i], point_x[j], rmj[j]);
			udm_sig_two[i] += Surface_gam[j] * deltaSn[j] * Umj[j];
			vdm_sig_two[i] += Surface_gam[j] * deltaSn[j] * Vmj[j];
		}
	//	printf("%0.4f, %0.4f\n", udm_sig_two[i], vdm_sig_two[i]);
	}
//********************************************************************************************************

//*************************************||   prepare velocity for Vortex convection  ||*************************************
	for(i=0;i<vortexNum;i++)									//Add two item
	{
		udm[i] = udm_sig_one[i] + udm_sig_two[i] + Uinf;
		vdm[i] = vdm_sig_one[i] + vdm_sig_two[i] + Vinf;
		fprintf(velocity, "%0.6f  %0.6f\n", udm[i], vdm[i]);
	//	printf("%0.6f  %0.6f\n", udm[i], vdm[i]);
	}
//********************************************************************************************************

//*************************************||             Vortex convection             ||*************************************		
	for(i=0;i<vortexNum;i++)									
	{
		point_x_next[i] = point_x[i] + (udm[i]) * deltaT;
		point_y_next[i] = point_y[i] + (vdm[i]) * deltaT;
	
	}
	for(i=0;i<vortexNum;i++)								
	{
		outputstreamL_x[i] = point_x_next[i];
		outputstreamL_y[i] = point_y_next[i];
		fprintf(driftPath, "%0.6f  %0.6f\n", outputstreamL_x[i], outputstreamL_y[i]);
	}









	return 0;
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

float cacu_K_Sm_Sn_Winf(float DS, float XM, float YM, float XN, float YN, float COS, float SIN)
{
	float temp_Kmn, temp_xm, temp_ym, temp_xn, temp_yn, temp_cos, temp_sin, temp_ds;
	temp_xm=XM; temp_ym=YM; temp_xn=XN; temp_yn=YN; temp_cos=COS; temp_sin=SIN; temp_ds=DS;
	if((temp_xm==temp_xn) && (temp_ym==temp_yn))
	{
		temp_Kmn=-0.5;
	}
	else
	{
		temp_Kmn=(temp_ds/(2*PI))*(((temp_ym-temp_yn)*temp_cos-(temp_xm-temp_xn)*temp_sin)/((float)pow((temp_xm-temp_xn),2)+(float)pow((temp_ym-temp_yn),2)));
	}
	return temp_Kmn;
}

float cacu_rhsm_Winf(float U, float COS, float V, float SIN)
{
	float temp_rhs, temp_u, temp_cos, temp_v, temp_sin;
	temp_u=U; temp_cos=COS; temp_v=V; temp_sin=SIN;
	temp_rhs=((-1)*temp_u*temp_cos)+((-1)*temp_v*temp_sin);
	return temp_rhs;
}

float cacu_K_Sm_Sn(float DS, float XM, float YM, float XN, float YN, float COS, float SIN)
{
	float temp_Kmn, temp_xm, temp_ym, temp_xn, temp_yn, temp_cos, temp_sin, temp_ds;
	temp_xm=XM; temp_ym=YM; temp_xn=XN; temp_yn=YN; temp_cos=COS; temp_sin=SIN; temp_ds=DS;
	if((temp_xm==temp_xn) && (temp_ym==temp_yn))
	{
		temp_Kmn=-0.5;
	}
	else
	{
		temp_Kmn=(temp_ds/(2*PI))*(((temp_ym-temp_yn)*temp_cos-(temp_xm-temp_xn)*temp_sin) / ((float)pow((temp_xm-temp_xn),2)+(float)pow((temp_ym-temp_yn),2)));
	}
	
	return temp_Kmn;
}

float cacu_rhsm(float U, float COS, float V, float SIN, float GAM)
{
	float temp_rhsm, temp_u, temp_cos, temp_v, temp_sin, temp_gam;
	temp_u=U; temp_cos=COS; temp_v=V; temp_sin=SIN; temp_gam=GAM; temp_gam=GAM;
	temp_rhsm=temp_gam*((temp_u*temp_cos)+(temp_v*temp_sin)+1.0);
	return temp_rhsm;
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

float cacu_Umn(float Xm, float Ym, float Xn, float Yn)
{
	float temp_U, tempXm, tempYm, tempXn, tempYn;
	tempXm=Xm; tempYm=Ym; tempXn=Xn; tempYn=Yn;
	temp_U=(1/(2*PI))*( (tempYm-tempYn) / ( (float)pow((tempXm-tempXn),2) + (float)pow((tempYm-tempYn),2) ));
	return temp_U;
}

float cacu_Vmn(float Xm, float Ym, float Xn, float Yn)
{
	float temp_V, tempXm, tempYm, tempXn, tempYn;
	tempXm=Xm; tempYm=Ym; tempXn=Xn; tempYn=Yn;
	temp_V=(-1/(2*PI))*((tempXm-tempXn) / ( (float)pow((tempXm-tempXn),2) + (float)pow((tempYm-tempYn),2) ));
	return temp_V;
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

void matrix_inversion_using_elementary_operation_Winf( float original[N][N], float inverse[N][N], int dim )
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
