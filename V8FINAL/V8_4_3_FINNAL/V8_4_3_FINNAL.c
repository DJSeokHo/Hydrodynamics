#include <stdio.h>
#include <math.h>

#define M 31				//number of surfaces
#define N 200				//Array size
#define PI 3.1415926535		//PI
#define r 0.6				//R of point J
#define a 0.5				//R of cycle
#define GAM 0.05				//GAM
#define Winf 1.0			//Winf
#define flow_angle 0.0		//flow angle

void draw(void);			//mian function

float init_x(float, float);                     //initialization  X of the circle
float init_y(float, float);                     //initialization  Y of the circle
float init_cp(float, float);                    //initialization  Control points of the circle
float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);						//trans radian to degree

float cacu_K_Sm_Sn_Winf(float, float, float, float, float, float, float);	//cacu K(Sm, Sn) for Winf
float cacu_rhsm_Winf(float, float, float, float);					//cacu rhsm Winf only

float cacu_K_Sm_Sn(float, float, float, float, float, float, float);	//cacu K(Sm, Sn) for Winf
float cacu_rhsm(float, float, float, float);					//cacu rhsm Winf only

void matrix_inversion_using_elementary_operation_Winf( float original[N][N], float inverse[N][N], int dim );		//Output GAM(s)
void matrix_inversion_using_elementary_operation( float original[N][N], float inverse[N][N], int dim );		//Output GAM(s)

float cacu_rm(float, float, float, float);      //Get Rm
float cacu_Um(float, float, float);             //Get Um
float cacu_Vm(float, float, float);             //Get Vm

float init_tempXY(float, float, float);       //Get xn and yn for each edge
int cacu_nsubs(float,float);					//Get number of sub-elements
float cacu_rn(float, float, float, float);		//Get Rn
float cacu_New_Um(float, int);					//Get New Um
float cacu_New_Vm(float, int);					//Get New Um

main()
{
	draw();
	return 0;
}

void draw(void)
{
	int i, j, time, nsubs[N]={0,};				//for Loop
	FILE *circle=fopen("circle.dat","wt");			//for geometry
	FILE *pivotalP=fopen("pivotalP.dat","wt");		//for geometry
	FILE *driftPath=fopen("driftPath.dat","wt");	//for driftPath

	float cycx[N]={0,}, cycy[N]={0,}, xmid[N]={0,}, ymid[N]={0,}, deltaSn[N]={0,}, angle, slope[N]={0,}, cosbetam[N]={0,}, sinbetam[N]={0,};		//init the cycle
	
	float Uinf, Vinf;	
	float Coeff[N]={0,};																//for Martensen
	float K_Sm_Sn_Winf[N][N]={0,}, inv_Kmn_Winf[N][N]={0,}, Surface_gam_Winf[N]={0,}, rhsM_Winf[N]={0,};	//for Martensen
	float K_Sm_Sn[N][N]={0,}, inv_Kmn[N][N]={0,}, Surface_gam[N]={0,}, rhsM[N]={0,};	//for Martensen
	float Surface_gam_add[N]={0,};

	float xja, yja, xjb, yjb;
	float rmj[N]={0,}, Umj[N]={0,}, Vmj[N]={0,};

	float ub, vb, ub_sig, vb_sig, deltaT;

	float tempX[N]={0,}, tempY[N]={0,}, tempsubs[N]={0,}, xn[N]={0,}, yn[N]={0,};

	float rnj[N]={0,};
	float sigmaU[N]={0,}, sigmaV[N]={0,};
	float NewUmj[N]={0,}, NewVmj[N]={0,};

//*************************************|| Input profile and other data preparation  ||*************************************
	
	Uinf = Winf*(float)cos(ang_to_rad(flow_angle));				//init Uniform flow
	Vinf = Winf*(float)sin(ang_to_rad(flow_angle));				//init Uniform flow

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

	xja = -1.2094;  yja = 0.0541;				//init point GAM J
	fprintf(driftPath, "%0.4f %0.4f\n", xja, yja);

	for(time=0;time<100;time++)
	{
//*************************************||           Caculate the Umj & Vmj          ||*************************************

		for(i=0;i<M-1;i++)
		{
			rmj[i] = cacu_rm(xmid[i], xja, ymid[i], yja);
			Umj[i] = cacu_Um(ymid[i], yja, rmj[i]);
			Vmj[i] = cacu_Vm(xmid[i], xja, rmj[i]);
		}

//*************************************||       Caculate sub-elements Umj & Vmj     ||*************************************
	
		for(i=0;i<M-1;i++)
		{
			nsubs[i] = cacu_nsubs(rmj[i], deltaSn[i]);
			tempsubs[i] = (1/(float)(nsubs[i]));	
			sigmaU[i] = 0.0; sigmaV[i] = 0.0;
			for(j=0; j<nsubs[i]; j++)
			{	
				if(nsubs[i]==1)
				{
					tempX[j] = cycx[i];
					tempY[j] = cycy[i];
					tempX[j+1] = cycx[i+1];
					tempY[j+1] = cycy[i+1];
				}
				else if( (nsubs[i]!=1) && (j!=(nsubs[i]-1)) )
				{
					tempX[0] = cycx[i];
					tempY[0] = cycy[i];		
					tempX[j+1] = init_tempXY(cycx[i+1], cycx[i], tempsubs[i]);
					tempY[j+1] = init_tempXY(cycy[i+1], cycy[i], tempsubs[i]);
					tempsubs[i] += (1/(float)(nsubs[i]));
					tempX[j+2] = cycx[i+1];
					tempY[j+2] = cycy[i+1];
				}
				else
					break;
			}
			for(j=0; j<nsubs[i]; j++)
			{
				xn[j] = init_cp(tempX[j+1], tempX[j]);
				yn[j] = init_cp(tempY[j+1], tempY[j]);
				rnj[j] = cacu_rn(xn[j], xja, yn[j], yja);
				sigmaU[i] += ((yn[j]-yja)/(float)pow(rnj[j],2));		
				sigmaV[i] += ((xn[j]-xja)/(float)pow(rnj[j],2));
			}
			NewUmj[i] = cacu_New_Um(sigmaU[i], nsubs[i]);
			NewVmj[i] = cacu_New_Vm(sigmaV[i], nsubs[i]);
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

//*************************************||        Martensen analysis Umj & Vmj       ||*************************************

		for(i=0;i<M-1;i++)																//Martensen analysis
		{
			for(j=0;j<M-1;j++)
			{
				K_Sm_Sn[i][j] = cacu_K_Sm_Sn(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);	//(XM, YM, XN, YN, COS, SIN)
			}
			rhsM[i] = cacu_rhsm(NewUmj[i], cosbetam[i], NewVmj[i], sinbetam[i]);	//(U, cos, V, sin)
		}

//*************************************||          Matrix inversion Winf only       ||*************************************

		matrix_inversion_using_elementary_operation_Winf(K_Sm_Sn_Winf,inv_Kmn_Winf,M-1);
		for(i=0;i<M-1;i++)
		{
			Surface_gam_Winf[i]=0.0;
			for(j=0;j<M-1;j++)
			{
				Surface_gam_Winf[i]+=inv_Kmn_Winf[i][j]*rhsM_Winf[j];				//Get surface GAM
			}
		}

//*************************************||          Matrix inversion Umj & Vmj       ||*************************************

		matrix_inversion_using_elementary_operation(K_Sm_Sn,inv_Kmn,M-1);
		for(i=0;i<M-1;i++)
		{
			Surface_gam[i]=0.0;
			for(j=0;j<M-1;j++)
			{
				Surface_gam[i]+=inv_Kmn[i][j]*rhsM[j];				//Get surface GAM
			}
		}

//*************************************||      Speed addition from Uinf & GAM J     ||*************************************

		for(i=0;i<M-1;i++)
		{
			Surface_gam_add[i] = Surface_gam[i] + Surface_gam_Winf[i];
		}

//*************************************||         Star Converction drift path       ||*************************************
	
		deltaT = 0.025;
		ub_sig = 0.0; vb_sig = 0.0;		
		for(i=0;i<M-1;i++)
		{
			ub_sig += Surface_gam_add[i] * deltaSn[i] * NewUmj[i];
			vb_sig += Surface_gam_add[i] * deltaSn[i] * NewVmj[i];
		}
		ub = (-1) * ub_sig + Uinf;
		vb = (-1) * vb_sig + Vinf;

		xjb = xja + ub * deltaT;
		yjb = yja + vb * deltaT;

		fprintf(driftPath, "%0.4f %0.4f\n", xjb, yjb);

		xja = xjb; yja = yjb;
	}
}



float init_x(float cta, float R)
{
	float temp_x,temp_r;
	temp_r=R;
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
	temp_Kmn=0; temp_xm=0; temp_ym=0; temp_xn=0; temp_yn=0; temp_cos=0; temp_sin=0; temp_coe=0;
}

float cacu_rhsm(float U, float COS, float V, float SIN)
{
	float temp_rhsm, temp_u, temp_cos, temp_v, temp_sin;
	temp_u=U; temp_cos=COS; temp_v=V; temp_sin=SIN;
	temp_rhsm=(-1)*GAM*((temp_u*temp_cos)+(temp_v*temp_sin));
	return temp_rhsm;
	temp_rhsm=0; temp_u=0; temp_cos=0; temp_v=0; temp_sin=0;
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
	float temp_rhsm, temp_u, temp_cos, temp_v, temp_sin;
	temp_u=U; temp_cos=COS; temp_v=V; temp_sin=SIN;
	temp_rhsm=(-1)*((temp_u*temp_cos)+(temp_v*temp_sin));
	return temp_rhsm;
	temp_rhsm=0; temp_u=0; temp_cos=0; temp_v=0; temp_sin=0;
}

float cacu_rm(float XM, float X_, float YM, float Y_)
{
	float temp_rm_, temp_xm_, temp_x_, temp_ym_, temp_y_;
	temp_xm_=XM; temp_x_=X_; temp_ym_=YM; temp_y_=Y_;
	temp_rm_ = (float)sqrt( (float)pow((temp_xm_-temp_x_),2) + (float)pow((temp_ym_-temp_y_),2) );
	return temp_rm_;
	temp_xm_=0; temp_x_=0;temp_ym_=0;temp_y_=0;temp_rm_=0;
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
	temp_ym_=0; temp_y_=0; temp_rm_=0; temp_Um_=0;
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
	temp_xm_=0; temp_x_=0; temp_rm_=0; temp_Vm_=0;
}

int cacu_nsubs(float RMJ, float DELTASM)
{
	float temp_nsubs_, temp_rmj_, temp_deltasm_;
	int n;
	temp_rmj_=RMJ; temp_deltasm_=DELTASM;
	temp_nsubs_=1+((2*temp_deltasm_)/temp_rmj_);
	
	n=(int)temp_nsubs_;
	return n;	
}

float init_tempXY(float V_1, float V, float SUBS)
{
	float temp_xy, temp_subs_, temp_v1_, temp_v_;
	temp_subs_ = SUBS; temp_v1_ = V_1; temp_v_ = V;
	if(temp_subs_ != 1)
	{
		temp_xy = temp_subs_ * (temp_v1_ - temp_v_) + temp_v_;
	}
	return temp_xy;
}

float cacu_rn(float XN, float XJ, float YN, float YJ)
{
	float temp_rn_, temp_xn_, temp_xj_, temp_yn_, temp_yj_;
	temp_xn_=XN; temp_xj_=XJ; temp_yn_=YN; temp_yj_=YJ;
	temp_rn_ = (float)sqrt( (float)pow((temp_xn_-temp_xj_),2) + (float)pow((temp_yn_-temp_yj_),2) );
	return temp_rn_;
}

float cacu_New_Um(float PLSTOTAL, int NSUB_)
{
	float temp_new_Um_, temp_tot_;
	int temp_nsubs_;
	temp_tot_=PLSTOTAL; temp_nsubs_=NSUB_;
	temp_new_Um_ = (1/(2*PI*temp_nsubs_))*PLSTOTAL;
	return temp_new_Um_;
}

float cacu_New_Vm(float PLSTOTAL, int NSUB_)
{
	float temp_new_Vm_, temp_tot_;
	int temp_nsubs_;
	temp_tot_=PLSTOTAL; temp_nsubs_=NSUB_;
	temp_new_Vm_ = (-1/(2*PI*temp_nsubs_))*PLSTOTAL;
	return temp_new_Vm_;
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