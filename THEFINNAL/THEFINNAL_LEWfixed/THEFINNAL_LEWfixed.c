#include <stdio.h>
#include <math.h>

#define M 9					//number of surfaces
#define N 20				//Array size
#define PI 3.1415926535		//PI
#define a 0.5				//R of cycle
#define Winf 1.0			//Winf
#define flow_angle 0.0		//flow angle
#define iter 1			//number of iterations

void draw(void);			//mian function

float init_x(float, float);                     //initialization  X of the circle
float init_y(float, float);                     //initialization  Y of the circle
float init_cp(float, float);                    //initialization  Control points of the circle
float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);						//trans radian to degree

float cacu_epsl(float);		//set shed points

float cacu_K_Sm_Sn_Winf(float, float, float, float, float, float, float);	//cacu K(Sm, Sn) for Winf
float cacu_rhsm_Winf(float, float, float, float);					//cacu rhsm Winf only

float cacu_K_Sm_Sn(float, float, float, float, float, float, float);	//cacu K(Sm, Sn) for MJ
float cacu_rhsm(float, float, float, float, float);					//cacu rhsm for MJ 

void matrix_inversion_using_elementary_operation_Winf( float original[N][N], float inverse[N][N], int dim );		//Output GAM(s)
void matrix_inversion_using_elementary_operation( float original[N][N], float inverse[N][N], int dim );		//Output GAM(s)

float cacu_rm(float, float, float, float);      //Get Rmj
float cacu_Um(float, float, float);             //Get Umj
float cacu_Vm(float, float, float);             //Get Vmj

float cacu_Umn(float, float, float, float);		//Get Umn. Here m & n means convection shed points only ?
float cacu_Vmn(float, float, float, float);		//Get Vmn. Here m & n means convection shed points only ?

main()
{
	draw();
	return 0;
}

void draw(void)
{
	int i, j, k, time, VortexNum, VortexIndex;			//for Loop

	FILE *circle=fopen("circle.dat","wt");				//for geometry
	FILE *pivotalP=fopen("pivotalP.dat","wt");			//for geometry
	FILE *driftPath=fopen("driftPath.dat","wt");		//for driftPath A
	FILE *shed=fopen("shed.dat","wt");

	FILE *mix_A_Winf=fopen("martixA_Winf.dat","wt");
	FILE *mix_B_Winf=fopen("martixB_Winf.dat","wt");
	FILE *mix_X_Winf=fopen("martixX_Winf.dat","wt");

	FILE *mix_A=fopen("martixA.dat","wt");
	FILE *mix_B=fopen("martixB.dat","wt");
	FILE *mix_X=fopen("martixX.dat","wt");

	FILE *velocity=fopen("velocity.dat","wt");
	FILE *vortexGAM=fopen("vortexGAM.dat","wt");

	float cycx[N]={0,}, cycy[N]={0,}, xmid[N]={0,}, ymid[N]={0,}, deltaSn[N]={0,}, angle, cosbetam[N]={0,}, sinbetam[N]={0,};		//init the cycle
	float x[N], y[N], tx[N], ty[N], nx[N], ny[N], offset;	//for shedding vortex

	float Uinf, Vinf;			//Uinf and Vinf
	float Coeff[N]={0,};																//for Martensen
	float K_Sm_Sn_Winf[N][N]={0,}, inv_Kmn_Winf[N][N]={0,}, Surface_gam_Winf[N]={0,}, rhsM_Winf[N]={0,};	//for Martensen
	float K_Sm_Sn[N][N]={0,}, inv_Kmn[N][N]={0,}, Surface_gam[N]={0,}, rhsM[N]={0,};	//for Martensen

	float rmj[N][N]={0,}, Umj[N][N]={0,}, Vmj[N][N]={0,};	//for induced-velocity
	float Umn[N]={0,}, Vmn[N]={0,};		//for self-induced-velocity

	float delta_GAM[N], deltaT;			//all shed points strength and delta T

	float shed_x[N]={0,}, shed_y[N]={0,};	//all shed vortices									

	float point_x[N]={0,}, point_y[N]={0,};		//all shed vortices								

	float point_x_next[N]={0,}, point_y_next[N]={0,};							//for drift path

	float sigma_rhsm_for_mj[N][N]={0,}, sigma_rhsm_mj[N]={0,};					//for K function's second sigma item

	float udm[N]={0,}, vdm[N]={0,};
	float udm_sig_one[N]={0,}, udm_sig_two[N]={0,}, vdm_sig_one[N]={0,}, vdm_sig_two[N]={0,};

	float outputstreamL_x[N]={0,}, outputstreamL_y[N]={0,};
	float check_one=0.0, check_two=0.0;									//check condtion

//*************************************|| Input profile and other data preparation  ||*************************************
	
	Uinf = Winf*(float)cos(ang_to_rad(flow_angle));				//init Uniform flow
	Vinf = Winf*(float)sin(ang_to_rad(flow_angle));				//init Uniform flow
	angle = 0.0;												//init angle of circle

	for(i = 0; i <= M-1; i++)														//////////////
	{																				//          //
		cycx[i] = init_x(ang_to_rad(angle),a);										//init cycle//
		cycy[i] = init_y(ang_to_rad(angle),a);										//          //
		fprintf(circle, "%0.6f  %0.6f\n", cycx[i], cycy[i]);						//          //
		angle += 360.0 / (M-1);														//////////////
	}																				

	for(i = 0; i < M - 1; i++)														//////////////
	{																				//          //
		xmid[i] = init_cp(cycx[i+1],cycx[i]);										// init mid //
		ymid[i] = init_cp(cycy[i+1],cycy[i]);										//   point  //
		fprintf(pivotalP, "%0.6f  %0.6f\n", xmid[i], ymid[i]);						//          //
		x[i] = cycx[i+1] - cycx[i];
		y[i] = cycy[i+1] - cycy[i];
	}																																						//////////////

	for(i = 0; i < M - 1; i++)														//element's deltaS
	{																				//Coefficient of Martensen
		deltaSn[i] = (float)sqrt((float)pow((cycx[i+1] - cycx[i]),2) + (float)pow((cycy[i+1] - cycy[i]),2));
		Coeff[i] = deltaSn[i]/(2*PI); 
		tx[i] = x[i] / deltaSn[i];
		ty[i] = y[i] / deltaSn[i];
		nx[i] = ty[i];
		ny[i] = (-1) * tx[i];
	}
	
	for(i = 0; i < M - 1; i++)														//////////////
	{																				//   init   //										// sinbetam //
		cosbetam[i] = tx[i];
		sinbetam[i] = ty[i];
	}	

//*************************************||        Martensen analysis Winf only       ||*************************************

	fprintf(mix_A_Winf, "Winf only, U = 1.0:\n\n");
	fprintf(mix_B_Winf, "Winf only, U = 1.0:\n\n");
	fprintf(mix_X_Winf, "Winf only, U = 1.0:\n\n");

	for(i=0;i<M-1;i++)																//Martensen analysis
	{
		for(j=0;j<M-1;j++)
		{
			K_Sm_Sn_Winf[i][j] = cacu_K_Sm_Sn_Winf(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[j],sinbetam[j]) + deltaSn[j];  // ==>deltaSn 더할 것  
			fprintf(mix_A_Winf, "%0.6f ", K_Sm_Sn_Winf[i][j]);
		}
		fprintf(mix_A_Winf, "\n");
		rhsM_Winf[i] = cacu_rhsm_Winf(Uinf, cosbetam[i], Vinf, sinbetam[i]) ;	//(U, cos, V, sin)
		fprintf(mix_B_Winf, "%0.6f\n", rhsM_Winf[i]);
	}

//*************************************||          Matrix inversion Winf only       ||*************************************

	matrix_inversion_using_elementary_operation_Winf(K_Sm_Sn_Winf,inv_Kmn_Winf,M-1);

     printf(" Gamma \n" ) ;
	for(i=0;i<M-1;i++)
	{
		Surface_gam_Winf[i]=0.0;
		for(j=0;j<M-1;j++)
		{
			Surface_gam_Winf[i]+=inv_Kmn_Winf[i][j]*rhsM_Winf[j];				//Get surface GAM
		}
		fprintf(mix_X_Winf, "%0.6f\n", Surface_gam_Winf[i]);
 //     printf( "%0.6f\n", Surface_gam_Winf[i]);
	}
// 여기 까지 OK
//*************************************||            Shed all free vortex           ||*************************************

	deltaT = 0.383;
	offset = 0.5;

	fprintf(vortexGAM, "t=1:\n");

	for(i=0;i<M-1;i++)
	{
		shed_x[i] = xmid[i] + nx[i] * offset * deltaSn[i];
		shed_y[i] = ymid[i] + ny[i] * offset * deltaSn[i];
		delta_GAM[i] = Surface_gam_Winf[i] * deltaSn[i];
		fprintf(shed, "%0.6f %0.6f %0.6f \n", shed_x[i], shed_y[i]);
	}

	for(time = 0; time < iter; time++)
	{
		VortexNum = (time + 1) * (M-1);			 //Vortex number follow time steps 
														 //Two vortices crate for per step
		for(i=0;i<M-1;i++)
		{
			point_x[time*(M-1)+i] = shed_x[i];
			point_y[time*(M-1)+i] = shed_y[i];
		}

		if(time!=0)
		{
			fprintf(vortexGAM, "\nt=%d:\n", time+1);
			for(i=0;i<M-1;i++)
			{
				delta_GAM[time*(M-1)+i] = Surface_gam[i] * deltaSn[i];
				fprintf(vortexGAM, "%0.6f\n", delta_GAM[time*(M-1)+i]);
			}
		}

//*************************************||           Caculate the Umj & Vmj          ||*************************************

		for(i=0;i<VortexNum;i++)			//i: number of vortices
		{
			for(j=0;j<M-1;j++)				//j: number of surfaces
			{
				rmj[i][j] = cacu_rm(xmid[j], point_x[i], ymid[j], point_y[i]);
				Umj[i][j] = cacu_Um(ymid[j], point_y[i], rmj[i][j]);
				Vmj[i][j] = cacu_Vm(xmid[j], point_x[i], rmj[i][j]);
				sigma_rhsm_for_mj[i][j] = cacu_rhsm(Umj[i][j], cosbetam[j], Vmj[i][j], sinbetam[j], delta_GAM[i]);
			}
		}

//*************************************|| Martensen analysis Umj & Vmj of shed point||*************************************
   
		for(i=0;i<M-1;i++)					//i: number of surfaces
		{
			sigma_rhsm_mj[i] = 0.0;
			for(j=0;j<VortexNum;j++)		//j: number of vortices
			{
				sigma_rhsm_mj[i] += sigma_rhsm_for_mj[j][i];
			}

		//	printf(" %0.6f \n", sigma_rhsm_mj[i]) ;   // 여기가 안 맞는다 위의 계산 과정 검토 해라. 식 (8.32) 로 계산 해라 
		}

		fprintf(mix_A, "\nt=%d:\n", time+1);
		fprintf(mix_B, "\nt=%d:\n", time+1);
		fprintf(mix_X, "\nt=%d:\n", time+1);

		for(i=0;i<M-1;i++)
		{
			rhsM[i] = rhsM_Winf[i] - sigma_rhsm_mj[i];
			fprintf(mix_B, "%0.6f\n", rhsM[i]);
			printf("%0.6f\n", rhsM[i]);
		}
		fprintf(mix_B, "\n");
		
		for(i=0;i<M-1;i++)																//Martensen analysis
		{
			for(j=0;j<M-1;j++)
			{
			//	K_Sm_Sn[i][j] = cacu_K_Sm_Sn(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);	//(XM, YM, XN, YN, COS, SIN)
				K_Sm_Sn[i][j] = cacu_K_Sm_Sn(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i])+deltaSn[j];	//(XM, YM, XN, YN, COS, SIN)
				fprintf(mix_A, "%0.6f  ", K_Sm_Sn[i][j]);
			}
			fprintf(mix_A, "\n");
		}
		fprintf(mix_A, "\n");

		matrix_inversion_using_elementary_operation(K_Sm_Sn,inv_Kmn,M-1);

		for(i=0;i<M-1;i++)
		{
			Surface_gam[i]=0.0;
			for(j=0;j<M-1;j++)
			{
				Surface_gam[i]+=inv_Kmn[i][j]*rhsM[j];				//Get surface GAM
			}
			fprintf(mix_X, "%0.6f\n", Surface_gam[i]);
		}
		fprintf(mix_X, "\n");

//*************************************||     Cacu the convection of Umn & Vmn      ||*************************************

		for(i=0;i<VortexNum;i++)
		{
			udm_sig_one[i] = 0.0; 
			vdm_sig_one[i] = 0.0;

			for(j=0;j<VortexNum;j++)
			{
				if(i!=j)
				{
					udm_sig_one[i] += delta_GAM[i] * cacu_Umn(point_x[i], point_y[i], point_x[j], point_y[j]);
					vdm_sig_one[i] += delta_GAM[i] * cacu_Vmn(point_x[i], point_y[i], point_x[j], point_y[j]);
				}
			}
		}

		for(i=0;i<VortexNum;i++)
		{
			udm_sig_two[i] = 0.0;
			vdm_sig_two[i] = 0.0;
		
			for(j=0;j<M-1;j++)									//init two item
			{
				udm_sig_two[i] += Surface_gam[j] * deltaSn[j] * Umj[i][j];
				vdm_sig_two[i] += Surface_gam[j] * deltaSn[j] * Vmj[i][j];
			}
		}

//*************************************||   prepare velocity for Vortex convection  ||*************************************
		
		fprintf(velocity, "\nt=%d:\n", time+1);

		for(i=0;i<VortexNum;i++)									//Add two item
		{
			if(time<0)
			{
				Vinf = 0.0;
			}
			else
			{
				Vinf = 0.0;
			}

			udm[i] = udm_sig_one[i] + udm_sig_two[i] + Uinf;
			vdm[i] = vdm_sig_one[i] + vdm_sig_two[i] + Vinf;
			fprintf(velocity, "%0.6f  %0.6f\n", udm[i], vdm[i]);
		}

//*************************************||             Vortex convection             ||*************************************
		
		for(i=0;i<VortexNum;i++)									
		{
			point_x_next[i] = point_x[i] + udm[i] * deltaT;
			point_y_next[i] = point_y[i] + vdm[i] * deltaT;
		}

//*************************************||             Output streamline             ||*************************************

		if(time == iter - 1)
		{
			for(i=0;i<VortexNum;i++)								
			{
				outputstreamL_x[i] = point_x_next[i];
				outputstreamL_y[i] = point_y_next[i];
			}
		}

		for(i=0;i<VortexNum;i++)									
		{
			point_x[i] = point_x_next[i];
			point_y[i] = point_y_next[i];
		}
		
		fprintf(driftPath, "\nt=%d:\n", time+1);

		if(time == iter - 1)
		{
			for(i=0;i<VortexNum;i++)
			{	
				fprintf(driftPath, "%0.6f  %0.6f\n", outputstreamL_x[i], outputstreamL_y[i]);
			}
		}

//*************************************||                Advance deltaT             ||*************************************

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

float cacu_epsl(float DELTASA)
{
	float temp_epsl_, temp_deltaSa_;
	temp_deltaSa_=DELTASA;
	temp_epsl_=0.5 * temp_deltaSa_;
	return temp_epsl_;
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

float cacu_rhsm(float U, float COS, float V, float SIN, float GAM)
{
	float temp_rhsm, temp_u, temp_cos, temp_v, temp_sin, temp_gam;
	temp_u=U; temp_cos=COS; temp_v=V; temp_sin=SIN; temp_gam=GAM; temp_gam=GAM;
//	temp_rhsm=temp_gam*((temp_u*temp_cos)+(temp_v*temp_sin));
	temp_rhsm=temp_gam*(1.0+(temp_u*temp_cos)+(temp_v*temp_sin));
	return temp_rhsm;
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
}

float cacu_rhsm_Winf(float U, float COS, float V, float SIN)
{
	float temp_rhsm, temp_u, temp_cos, temp_v, temp_sin;
	temp_u=U; temp_cos=COS; temp_v=V; temp_sin=SIN;
	temp_rhsm=(-1)*((temp_u*temp_cos)+(temp_v*temp_sin));
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