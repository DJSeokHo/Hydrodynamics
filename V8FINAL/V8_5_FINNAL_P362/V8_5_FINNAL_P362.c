#include <stdio.h>
#include <math.h>

#define M 9				//number of surfaces
#define N 500			//Array size
#define PI 3.1415926535		//PI
#define a 0.5				//R of cycle
#define Winf 1.0			//Winf
#define flow_angle 0.0		//flow angle
#define iter 1				//number of iterations

void draw(void);			//mian function

float init_x(float, float);                     //initialization  X of the circle
float init_y(float, float);                     //initialization  Y of the circle
float init_cp(float, float);                    //initialization  Control points of the circle
float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);						//trans radian to degree

float cacu_deltaGAM(float gam, float deltaT);	//set shed points
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
	int i, j, k, time, VortexNum, VortexIndex;				//for Loop
	FILE *circle=fopen("circle.dat","wt");			//for geometry
	FILE *pivotalP=fopen("pivotalP.dat","wt");		//for geometry
	FILE *driftPathA=fopen("driftPathA.dat","wt");	//for driftPath A
	FILE *driftPathB=fopen("driftPathB.dat","wt");	//for driftPath B
	FILE *shedA=fopen("shedA.dat","wt");				//for shed points A
	FILE *shedB=fopen("shedB.dat","wt");				//for shed points B
	FILE *mtixA=fopen("matrix_A.dat","wt");
	FILE *mtixB_before=fopen("matrix_B_before.dat","wt");
	FILE *mtixB_after=fopen("matrix_B_after.dat","wt");
	FILE *mtixX=fopen("matrix_X.dat","wt");

	float cycx[N]={0,}, cycy[N]={0,}, xmid[N]={0,}, ymid[N]={0,}, deltaSn[N]={0,}, angle, slope[N]={0,}, cosbetam[N]={0,}, sinbetam[N]={0,};		//init the cycle
	
	float Uinf, Vinf;	
	float Coeff[N]={0,};																//for Martensen
	float K_Sm_Sn_Winf[N][N]={0,}, inv_Kmn_Winf[N][N]={0,}, Surface_gam_Winf[N]={0,}, rhsM_Winf[N]={0,};	//for Martensen
	float K_Sm_Sn[N][N]={0,}, inv_Kmn[N][N]={0,}, Surface_gam[N]={0,}, rhsM[N]={0,};	//for Martensen
	float Surface_gam_add[N]={0,};														//for Martensen

	float rmj[N]={0,}, Umj[N]={0,}, Vmj[N]={0,};
	float Umn[N]={0,}, Vmn[N]={0,};

	float delta_GAM_A, delta_GAM_B, deltaT, epsl_length, ceta;			//two shed points
	float shed_A_x, shed_A_y;															//two shed points
	float shed_B_x, shed_B_y;															//two shed points
	float point_A_x[N]={0,}, point_A_y[N]={0,};											//two shed points
	float point_B_x[N]={0,}, point_B_y[N]={0,};											//two shed points
	float point_A_x_next[N]={0,}, point_A_y_next[N]={0,};												//for drift path
	float point_B_x_next[N]={0,}, point_B_y_next[N]={0,};								//for drift path
	
	float sigma_rhsm_for_mj[N][N]={0,}, sigma_rhsm_mj[N]={0,};					//for K function's second sigma item
	
	float udm[N]={0,}, vdm[N]={0,};
	float udm_sig_one[N]={0,}, udm_sig_two[N]={0,}, vdm_sig_one[N]={0,}, vdm_sig_two[N]={0,};

	float outputstreamL_A_x[N]={0,}, outputstreamL_A_y[N]={0,};
	float outputstreamL_B_x[N]={0,}, outputstreamL_B_y[N]={0,};

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
	}																				

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
		cosbetam[i] = ((cycx[i] - cycx[i+1]))/deltaSn[i];							//////////////
		sinbetam[i] = ((cycy[i] - cycy[i+1]))/deltaSn[i];							
	}	

//*************************************||        Martensen analysis Winf only       ||*************************************

	for(i=0;i<M-1;i++)																//Martensen analysis
	{
		for(j=0;j<M-1;j++)
		{
			K_Sm_Sn_Winf[i][j] = cacu_K_Sm_Sn_Winf(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);	//(XM, YM, XN, YN, COS, SIN)
		//	fprintf(mtixA, "%8.4f  ", K_Sm_Sn_Winf[i][j]);
		}
	//	fprintf(mtixA, "\n", K_Sm_Sn_Winf[i][j]);
		rhsM_Winf[i] = cacu_rhsm_Winf(Uinf, cosbetam[i], Vinf, sinbetam[i]);	//(U, cos, V, sin)
		fprintf(mtixB_before, "%8.4f\n", rhsM_Winf[i]);
	//	printf("%d:  rhsM_Winf is %0.4f\n", i+1, rhsM_Winf[i]);
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
	//	fprintf(mtixX, "%8.4f\n", Surface_gam_Winf[i]);
	//	printf("%d:  Surface_gam_Winf is %0.4f\n", i+1, Surface_gam_Winf[i]);
	}

	deltaT = 0.05;																	//T is using for make shed point
	
//*************************************||            Shed two free vortex           ||*************************************
	
	epsl_length = cacu_epsl(deltaSn[2]);											//init shed points
	ceta = 45;	
	shed_A_x=epsl_length*(float)cos(ang_to_rad(ceta))+xmid[2];
	shed_A_y=epsl_length*(float)sin(ang_to_rad(ceta))+ymid[2];

	fprintf(shedA, "%0.4f %0.4f\n", shed_A_x, shed_A_y);
	
	outputstreamL_A_x[iter+1] = shed_A_x;
	outputstreamL_A_y[iter+1] = shed_A_y;

	delta_GAM_A = cacu_deltaGAM(Surface_gam_Winf[3], deltaT);

	for(time = 0; time < iter; time++)
	{
		VortexNum = time + 1;
		point_A_x[time] = shed_A_x;
		point_A_y[time] = shed_A_y;

//*************************************||           Caculate the Umj & Vmj          ||*************************************
	
		for(i=0;i<VortexNum;i++)
		{
			for(j=0;j<M-1;j++)	//Here means only have one free vortex
			{
				rmj[j] = cacu_rm(xmid[j], point_A_x[i], ymid[j], point_A_y[i]);
				Umj[j] = cacu_Um(ymid[j], point_A_y[i], rmj[j]);
				Vmj[j] = cacu_Vm(xmid[j], point_A_x[i], rmj[j]);
			//	printf("%d:  rmj is %0.4f   U,V is %0.4f   %0.4f\n", j+1, rmj[j], Umj[j], Vmj[j]);
				sigma_rhsm_for_mj[i][j] = cacu_rhsm(Umj[j], cosbetam[j], Vmj[j], sinbetam[j], delta_GAM_A);
			//	printf("%d:  sigma rhsm or mj is %0.4f\n", i+1, sigma_rhsm_for_mj[i][j]);
			}
		}

//*************************************|| Martensen analysis Umj & Vmj of shed point||*************************************

		for(i=0;i<M-1;i++)
		{
			sigma_rhsm_mj[i] = 0.0;
			for(j=0;j<VortexNum;j++)
			{
				sigma_rhsm_mj[i] += sigma_rhsm_for_mj[j][i];
			}
		//	printf("%d:  sigma rhsm or mj is %0.4f\n", i+1, sigma_rhsm_mj[i]);
		}

		for(i=0;i<M-1;i++)
		{
			rhsM[i] = rhsM_Winf[i] - sigma_rhsm_mj[i];
		//	printf("%d:  rhsm is %0.4f\n", i+1, rhsM[i]);
		}
		
		for(i=0;i<M-1;i++)																//Martensen analysis
		{
			for(j=0;j<M-1;j++)
			{
				K_Sm_Sn[i][j] = cacu_K_Sm_Sn(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]) + deltaSn[i];	//(XM, YM, XN, YN, COS, SIN)
			//	K_Sm_Sn[i][j] = cacu_K_Sm_Sn(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);
			}
		}

		matrix_inversion_using_elementary_operation_Winf(K_Sm_Sn,inv_Kmn,M-1);
		for(i=0;i<M-1;i++)
		{
			Surface_gam[i]=0.0;
			for(j=0;j<M-1;j++)
			{
				Surface_gam[i]+=inv_Kmn[i][j]*rhsM[j];				//Get surface GAM
			}
		//	printf("%d:  Surface_gam is %0.4f\n", i+1, Surface_gam[i]);
		}

//*************************************||     Cacu the convection of Umn & Vmn      ||*************************************

		if(VortexNum>1)
		{
			for(i=0;i<VortexNum;i++)									//init two item
			{
				udm_sig_one[i] = 0.0; vdm_sig_one[i] = 0.0;
				for(j=0;j<VortexNum;j++)
				{
					if(i!=j)
					{
						udm_sig_one[i] += delta_GAM_A * cacu_Umn(point_A_x[i], point_A_y[i], point_A_x[j], point_A_y[j]);
						vdm_sig_one[i] += delta_GAM_A * cacu_Vmn(point_A_x[i], point_A_y[i], point_A_x[j], point_A_y[j]);
					}	
				}
			}
		}

		for(i=0;i<VortexNum;i++)									//init two item
		{
			udm_sig_two[i] = 0.0; vdm_sig_two[i] = 0.0;
			for(j=0;j<M-1;j++)
			{
				udm_sig_two[i] += Surface_gam[j] * deltaSn[j] * Umj[j];
				vdm_sig_two[i] += Surface_gam[j] * deltaSn[j] * Vmj[j];
			}
		}

//*************************************||   prepare velocity for Vortex convection  ||*************************************
		
		for(i=0;i<VortexNum;i++)									//init two item
		{
			udm[i] = udm_sig_one[i] + udm_sig_two[i];
			vdm[i] = vdm_sig_one[i] + vdm_sig_two[i];
		}

//*************************************||             Vortex convection             ||*************************************
		
		for(i=0;i<VortexNum;i++)									//init two item
		{
			point_A_x_next[i] = point_A_x[i] + (udm[i] + Uinf) * deltaT;
			point_A_y_next[i] = point_A_y[i] + (vdm[i] + Vinf) * deltaT;
		}
		
//*************************************||             Output streamline             ||*************************************

		if(time == iter - 1)
		{
			for(i=0;i<VortexNum;i++)									//init two item
			{
			//	fprintf(driftPathA, "%0.4f %0.4f\n", point_A_x_next[i], point_A_y_next[i]);
				outputstreamL_A_x[i+1] = point_A_x_next[i];
				outputstreamL_A_y[i+1] = point_A_y_next[i];
			}
		}

		for(i=0;i<VortexNum;i++)									//init two item
		{
			point_A_x[i] = point_A_x_next[i];
			point_A_y[i] = point_A_y_next[i];
		}
		
		if(time == iter - 1)
		{
		//	fprintf( driftPathA, "\nTitle='streamline'\n");
		//	fprintf( driftPathA, "Variables='x','y'\n");
		//	fprintf( driftPathA, "Zone I=%d\n\n", VortexNum);
			for(i=VortexNum+1;i>0;i--)
			{
				fprintf(driftPathA, "%0.4f %0.4f\n", outputstreamL_A_x[i], outputstreamL_A_y[i]);
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

float cacu_deltaGAM(float gam, float deltaT)
{
	float temp_deltaGAM_, temp_gam_, temp_deltaT_;
	temp_gam_=gam; temp_deltaT_=deltaT;
	temp_deltaGAM_=0.5 * (float)pow(temp_gam_,2)*temp_deltaT_;
	return temp_deltaGAM_;
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
	temp_rhsm=temp_gam * ( 1 + (temp_u*temp_cos)+(temp_v*temp_sin) );
//	temp_rhsm=temp_gam * ((temp_u*temp_cos)+(temp_v*temp_sin));
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