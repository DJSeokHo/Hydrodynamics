#include <stdio.h>
#include <math.h>

#define M 11				//number of surfaces
#define N 200				//Array size
#define PI 3.1415926535		//PI
#define a 0.5				//R of cycle
#define Winf 1.0			//Winf
#define flow_angle 0.0		//flow angle
#define iterA 20			//number of iterations
#define iterB 20			//number of iterations
#define order 1				//number of central difference order
#define FVN 2				//number of free vortex shed places

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

	float cycx[N]={0,}, cycy[N]={0,}, xmid[N]={0,}, ymid[N]={0,}, deltaSn[N]={0,}, angle, slope[N]={0,}, cosbetam[N]={0,}, sinbetam[N]={0,};		//init the cycle
	
	int shedIndex[N]={0,};

	float Uinf, Vinf;	
	float Coeff[N]={0,};																//for Martensen

	float K_Sm_Sn_Winf[N][N]={0,}, inv_Kmn_Winf[N][N]={0,}, Surface_gam_Winf[N]={0,}, rhsM_Winf[N]={0,};	//for Martensen

	float K_Sm_Sn_A[N][N]={0,}, inv_Kmn_A[N][N]={0,}, Surface_gam_A[N]={0,}, rhsM_A[N]={0,};	//for Martensen
	float Surface_gam_add_A[N]={0,};														//for Martensen

	float K_Sm_Sn_B[N][N]={0,}, inv_Kmn_B[N][N]={0,}, Surface_gam_B[N]={0,}, rhsM_B[N]={0,};	//for Martensen
	float Surface_gam_add_B[N]={0,};														//for Martensen

	float rmj_A[N]={0,}, Umj_A[N]={0,}, Vmj_A[N]={0,};
	float Umn_A[N]={0,}, Vmn_A[N]={0,};

	float rmj_B[N]={0,}, Umj_B[N]={0,}, Vmj_B[N]={0,};
	float Umn_B[N]={0,}, Vmn_B[N]={0,};

	float delta_GAM_A, delta_GAM_B, deltaT, epsl_length, ceta;			//two shed points
	float shed_A_x, shed_A_y;															//two shed points
	float shed_B_x, shed_B_y;															//two shed points
	float point_A_x[N]={0,}, point_A_y[N]={0,};											//two shed points
	float point_B_x[N]={0,}, point_B_y[N]={0,};											//two shed points
	float point_A_x_next[N]={0,}, point_A_y_next[N]={0,};								//for drift path
	float point_B_x_next[N]={0,}, point_B_y_next[N]={0,};								//for drift path
	
	float sigma_rhsm_for_mj_A[N][N]={0,}, sigma_rhsm_mj_A[N]={0,};					//for K function's second sigma item
	float sigma_rhsm_for_mj_B[N][N]={0,}, sigma_rhsm_mj_B[N]={0,};					//for K function's second sigma item
	
	float udm_A[N]={0,}, vdm_A[N]={0,};
	float udm_sig_one_A[N]={0,}, udm_sig_two_A[N]={0,}, vdm_sig_one_A[N]={0,}, vdm_sig_two_A[N]={0,};

	float udm_B[N]={0,}, vdm_B[N]={0,};
	float udm_sig_one_B[N]={0,}, udm_sig_two_B[N]={0,}, vdm_sig_one_B[N]={0,}, vdm_sig_two_B[N]={0,};

	float outputstreamL_A_x[N]={0,}, outputstreamL_A_y[N]={0,};
	float outputstreamL_B_x[N]={0,}, outputstreamL_B_y[N]={0,};

//*************************************|| Input profile and other data preparation  ||*************************************
	
	Uinf = Winf*(float)cos(ang_to_rad(flow_angle));				//init Uniform flow
	Vinf = Winf*(float)sin(ang_to_rad(flow_angle));				//init Uniform flow
	angle = 0.0;												//init angle of circle
	
	for(i = 0; i <= M-1; i++)	//M is number of point								//////////////
	{																				//          //
		cycx[i] = init_x(ang_to_rad(angle),a);										//init cycle//
		cycy[i] = init_y(ang_to_rad(angle),a);										//          //
		fprintf(circle, "%0.4f  %0.4f\n", cycx[i], cycy[i]);						//          //
		angle += 360.0 / (M-1);														//////////////
	}																				

	for(i = 0; i < M - 1; i++)	//M is number of point								//////////////
	{																				//          //
		xmid[i] = init_cp(cycx[i+1],cycx[i]);										// init mid //
		ymid[i] = init_cp(cycy[i+1],cycy[i]);										//   point  //
		fprintf(pivotalP, "%0.4f  %0.4f\n", xmid[i], ymid[i]);						//          //
	}																				//////////////		


	for(i = 0; i < M - 1; i++)	//M is number of point								//element's deltaS
	{																				//Coefficient of Martensen
		deltaSn[i] = (float)sqrt((float)pow((cycx[i+1] - cycx[i]),2) + (float)pow((cycy[i+1] - cycy[i]),2));
		Coeff[i] = deltaSn[i]/(2*PI);
	}

	for(i = 0; i < M - 1; i++)	//M is number of point								//////////////
	{																				//   init   //										// sinbetam //
		cosbetam[i] = ((cycx[i] - cycx[i+1]))/deltaSn[i];							//////////////
		sinbetam[i] = ((cycy[i] - cycy[i+1]))/deltaSn[i];
		if(cosbetam[i] > 0.999999)
		{
			shedIndex[0] = i;
		}
		if(cosbetam[i] < -0.99999)
		{
			shedIndex[1] = i;
		}
	//	printf("%0.8f  %0.8f\n", cosbetam[i], sinbetam[i]);
	}
//	printf("%d  %d\n", shedIndex[0], shedIndex[1]);

//*************************************||        Martensen analysis Winf only       ||*************************************

	for(i=0;i<M-1;i++)		//M is number of point									//Martensen analysis
	{
		for(j=0;j<M-1;j++)     //M is number of point	  
		{
			K_Sm_Sn_Winf[i][j] = cacu_K_Sm_Sn_Winf(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);	//(XM, YM, XN, YN, COS, SIN)
		}
		rhsM_Winf[i] = cacu_rhsm_Winf(Uinf, cosbetam[i], Vinf, sinbetam[i]);	//(U, cos, V, sin)
	//	printf("%d -- %0.4f\n", i+1, rhsM_Winf[i]);
	}

//*************************************||          Matrix inversion Winf only       ||*************************************

	matrix_inversion_using_elementary_operation_Winf(K_Sm_Sn_Winf,inv_Kmn_Winf,M-1);
	for(i=0;i<M-1;i++)		//M is number of point	
	{
		Surface_gam_Winf[i]=0.0;
		for(j=0;j<M-1;j++)		//M is number of point	
		{
			Surface_gam_Winf[i]+=inv_Kmn_Winf[i][j]*rhsM_Winf[j];				//Get surface GAM
		}
	//	printf("%d -- %0.4f\n", i+1, Surface_gam_Winf[i]);
	}

	deltaT = 0.05;																	//T is using for make shed point
	
//*************************************||            Shed two free vortex           ||*************************************
	
	epsl_length = 0.05;											//init shed points
	ceta = 45;	
	shed_A_x=epsl_length*(float)cos(ang_to_rad(ceta))+xmid[shedIndex[0]];
	shed_A_y=epsl_length*(float)sin(ang_to_rad(ceta))+ymid[shedIndex[0]];
	shed_B_x=epsl_length*(float)cos(ang_to_rad((-1)*ceta))+xmid[shedIndex[1]];
	shed_B_y=epsl_length*(float)sin(ang_to_rad((-1)*ceta))+ymid[shedIndex[1]];

	fprintf(shedA, "%0.4f %0.4f\n", shed_A_x, shed_A_y);
	fprintf(shedB, "%0.4f %0.4f\n", shed_B_x, shed_B_y);

	outputstreamL_A_x[iterA+1] = shed_A_x;
	outputstreamL_A_y[iterA+1] = shed_A_y;

	outputstreamL_B_x[iterB+1] = shed_B_x;
	outputstreamL_B_y[iterB+1] = shed_B_y;

	delta_GAM_A = cacu_deltaGAM(Surface_gam_Winf[shedIndex[0]+1], deltaT);
	delta_GAM_B = (-1) * cacu_deltaGAM(Surface_gam_Winf[shedIndex[1]-1], deltaT);

//	printf("%0.4f  %0.4f\n", delta_GAM_A, delta_GAM_B);

	for(time = 0; time < iterA; time++)		//time steps
	{
		VortexNum = time + 1;
		point_A_x[time] = shed_A_x;
		point_A_y[time] = shed_A_y;

//*************************************||           Caculate the Umj & Vmj          ||*************************************

		for(i=0;i<VortexNum;i++)
		{
			for(j=0;j<M-1;j++)	//Here means only have one free vortex
			{
				rmj_A[j] = cacu_rm(xmid[j], point_A_x[i], ymid[j], point_A_y[i]);
				Umj_A[j] = cacu_Um(ymid[j], point_A_y[i], rmj_A[j]);
				Vmj_A[j] = cacu_Vm(xmid[j], point_A_x[i], rmj_A[j]);
				sigma_rhsm_for_mj_A[i][j] = cacu_rhsm(Umj_A[j], cosbetam[j], Vmj_A[j], sinbetam[j], delta_GAM_A);
			}
		}

//*************************************|| Martensen analysis Umj & Vmj of shed point||*************************************

		for(i=0;i<M-1;i++)
		{
			sigma_rhsm_mj_A[i] = 0.0;
			for(j=0;j<VortexNum;j++)
			{
				sigma_rhsm_mj_A[i] += sigma_rhsm_for_mj_A[j][i];
			}
		}

		for(i=0;i<M-1;i++)
		{
			rhsM_A[i] = rhsM_Winf[i] - sigma_rhsm_mj_A[i];
		}
		
		for(i=0;i<M-1;i++)																//Martensen analysis
		{
			for(j=0;j<M-1;j++)
			{
				K_Sm_Sn_A[i][j] = cacu_K_Sm_Sn(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);	//(XM, YM, XN, YN, COS, SIN)
			}
		}

		matrix_inversion_using_elementary_operation(K_Sm_Sn_A,inv_Kmn_A,M-1);
		for(i=0;i<M-1;i++)
		{
			Surface_gam_A[i]=0.0;
			for(j=0;j<M-1;j++)
			{
				Surface_gam_A[i]+=inv_Kmn_A[i][j]*rhsM_A[j];				//Get surface GAM
			}
		}

//*************************************||     Cacu the convection of Umn & Vmn      ||*************************************

		if(VortexNum>1)
		{
			for(i=0;i<VortexNum;i++)									//init two item
			{
				udm_sig_one_A[i] = 0.0; vdm_sig_one_A[i] = 0.0;
				for(j=0;j<VortexNum;j++)
				{
					if(i!=j)
					{
						udm_sig_one_A[i] += delta_GAM_A * cacu_Umn(point_A_x[i], point_A_y[i], point_A_x[j], point_A_y[j]);
						vdm_sig_one_A[i] += delta_GAM_A * cacu_Vmn(point_A_x[i], point_A_y[i], point_A_x[j], point_A_y[j]);
					}	
				}
			}
		}

		for(i=0;i<VortexNum;i++)									//init two item
		{
			udm_sig_two_A[i] = 0.0; vdm_sig_two_A[i] = 0.0;
			for(j=0;j<M-1;j++)
			{
				udm_sig_two_A[i] += Surface_gam_A[j] * deltaSn[j] * Umj_A[j];
				vdm_sig_two_A[i] += Surface_gam_A[j] * deltaSn[j] * Vmj_A[j];
			}
		}

//*************************************||   prepare velocity for Vortex convection  ||*************************************
		
		for(i=0;i<VortexNum;i++)									//init two item
		{
			udm_A[i] = udm_sig_one_A[i] + udm_sig_two_A[i];
			vdm_A[i] = vdm_sig_one_A[i] + vdm_sig_two_A[i];
		}

//*************************************||             Vortex convection             ||*************************************
		
		for(i=0;i<VortexNum;i++)									//init two item
		{
			point_A_x_next[i] = point_A_x[i] + (udm_A[i] + Uinf) * deltaT;
			point_A_y_next[i] = point_A_y[i] + (vdm_A[i] + Vinf) * deltaT;
			tempCD_x[i] = point_A_x_next[i];
			tempCD_y[i] = point_A_y_next[i];
		} 

//*************************************||             Output streamline             ||*************************************

		if(time == iterA - 1)
		{
			for(i=0;i<VortexNum;i++)									//init two item
			{
				outputstreamL_A_x[i+1] = point_A_x_next[i];
				outputstreamL_A_y[i+1] = point_A_y_next[i];
			}
		}
		

		for(i=0;i<VortexNum;i++)									//init two item
		{
			point_A_x[i] = point_A_x_next[i];
			point_A_y[i] = point_A_y_next[i];
		}
		
		if(time == iterA - 1)
		{
			for(i=1;i<VortexNum+2;i++)
			{
				fprintf(driftPathA, "%0.4f %0.4f\n", outputstreamL_A_x[i], outputstreamL_A_y[i]);
			}
		}

//*************************************||                Advance deltaT             ||*************************************

	}






	for(time = 0; time < iterB; time++)		//time steps
	{
		VortexNum = time + 1;
		point_B_x[time] = shed_B_x;
		point_B_y[time] = shed_B_y;

//*************************************||           Caculate the Umj & Vmj          ||*************************************
	
		for(i=0;i<VortexNum;i++)
		{
			for(j=0;j<M-1;j++)	//Here means only have one free vortex
			{
				rmj_B[j] = cacu_rm(xmid[j], point_B_x[i], ymid[j], point_B_y[i]);
				Umj_B[j] = cacu_Um(ymid[j], point_B_y[i], rmj_B[j]);
				Vmj_B[j] = cacu_Vm(xmid[j], point_B_x[i], rmj_B[j]);
				sigma_rhsm_for_mj_B[i][j] = cacu_rhsm(Umj_B[j], cosbetam[j], Vmj_B[j], sinbetam[j], delta_GAM_B);
			}
		}

//*************************************|| Martensen analysis Umj & Vmj of shed point||*************************************

		for(i=0;i<M-1;i++)
		{
			sigma_rhsm_mj_B[i] = 0.0;
			for(j=0;j<VortexNum;j++)
			{
				sigma_rhsm_mj_B[i] += sigma_rhsm_for_mj_B[j][i];
			}
		}

		for(i=0;i<M-1;i++)
		{
			rhsM_B[i] = rhsM_Winf[i] - sigma_rhsm_mj_B[i];
		}
		
		for(i=0;i<M-1;i++)																//Martensen analysis
		{
			for(j=0;j<M-1;j++)
			{
				K_Sm_Sn_B[i][j] = cacu_K_Sm_Sn(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);	//(XM, YM, XN, YN, COS, SIN)
			}
		}

		matrix_inversion_using_elementary_operation(K_Sm_Sn_B,inv_Kmn_B,M-1);
		for(i=0;i<M-1;i++)
		{
			Surface_gam_B[i]=0.0;
			for(j=0;j<M-1;j++)
			{
				Surface_gam_B[i]+=inv_Kmn_B[i][j]*rhsM_B[j];				//Get surface GAM
			}
		}

//*************************************||     Cacu the convection of Umn & Vmn      ||*************************************

		if(VortexNum>1)
		{
			for(i=0;i<VortexNum;i++)									//init two item
			{
				udm_sig_one_B[i] = 0.0; vdm_sig_one_B[i] = 0.0;
				for(j=0;j<VortexNum;j++)
				{
					if(i!=j)
					{
						udm_sig_one_B[i] += delta_GAM_B * cacu_Umn(point_B_x[i], point_B_y[i], point_B_x[j], point_B_y[j]);
						vdm_sig_one_B[i] += delta_GAM_B * cacu_Vmn(point_B_x[i], point_B_y[i], point_B_x[j], point_B_y[j]);
					}	
				}
			}
		}

		for(i=0;i<VortexNum;i++)									//init two item
		{
			udm_sig_two_B[i] = 0.0; vdm_sig_two_B[i] = 0.0;
			for(j=0;j<M-1;j++)
			{
				udm_sig_two_B[i] += Surface_gam_B[j] * deltaSn[j] * Umj_B[j];
				vdm_sig_two_B[i] += Surface_gam_B[j] * deltaSn[j] * Vmj_B[j];
			}
		}

//*************************************||   prepare velocity for Vortex convection  ||*************************************
		
		for(i=0;i<VortexNum;i++)									//init two item
		{
			udm_B[i] = udm_sig_one_B[i] + udm_sig_two_B[i];
			vdm_B[i] = vdm_sig_one_B[i] + vdm_sig_two_B[i];
		}

//*************************************||             Vortex convection             ||*************************************
		
		for(i=0;i<VortexNum;i++)									//init two item
		{
			point_B_x_next[i] = point_B_x[i] + (udm_B[i] + Uinf) * deltaT;
			point_B_y_next[i] = point_B_y[i] + (vdm_B[i] + Vinf) * deltaT;
		}
		
//*************************************||             Output streamline             ||*************************************

		if(time == iterB - 1)
		{
			for(i=0;i<VortexNum;i++)									//init two item
			{
				outputstreamL_B_x[i+1] = point_B_x_next[i];
				outputstreamL_B_y[i+1] = point_B_y_next[i];
			}
		}

		for(i=0;i<VortexNum;i++)									//init two item
		{
			point_B_x[i] = point_B_x_next[i];
			point_B_y[i] = point_B_y_next[i];
		}
		
		if(time == iterB - 1)
		{
			for(i=0;i<VortexNum+2;i++)
			{
				fprintf(driftPathB, "%0.4f %0.4f\n", outputstreamL_B_x[i], outputstreamL_B_y[i]);
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

float cacu_deltaGAM(float gam, float deltaT)
{
	float temp_deltaGAM_, temp_gam_, temp_deltaT_;
	temp_gam_=gam; temp_deltaT_=deltaT;
	temp_deltaGAM_=0.5 * (float)pow(temp_gam_,2)*temp_deltaT_;
	return temp_deltaGAM_;
	temp_deltaGAM_=0; temp_gam_=0; temp_deltaT_=0;
}

float cacu_epsl(float DELTASA)
{
	float temp_epsl_, temp_deltaSa_;
	temp_deltaSa_=DELTASA;
	temp_epsl_=0.5 * temp_deltaSa_;
	return temp_epsl_;
	temp_epsl_=0; temp_deltaSa_=0;
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

float cacu_rhsm(float U, float COS, float V, float SIN, float GAM)
{
	float temp_rhsm, temp_u, temp_cos, temp_v, temp_sin, temp_gam;
	temp_u=U; temp_cos=COS; temp_v=V; temp_sin=SIN; temp_gam=GAM; temp_gam=GAM;
	temp_rhsm=temp_gam*((temp_u*temp_cos)+(temp_v*temp_sin));
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

float cacu_Umn(float Xm, float Ym, float Xn, float Yn)
{
	float temp_U, tempXm, tempYm, tempXn, tempYn;
	tempXm=Xm; tempYm=Ym; tempXn=Xn; tempYn=Yn;
	temp_U=(1/(2*PI))*( (tempYm-tempYn) / ( (float)pow((tempXm-tempXn),2) + (float)pow((tempYm-tempYn),2) ));
	return temp_U;
	temp_U=0; tempXm=0; tempYm=0; tempXn=0; tempYn=0;
}

float cacu_Vmn(float Xm, float Ym, float Xn, float Yn)
{
	float temp_V, tempXm, tempYm, tempXn, tempYn;
	tempXm=Xm; tempYm=Ym; tempXn=Xn; tempYn=Yn;
	temp_V=(-1/(2*PI))*((tempXm-tempXn) / ( (float)pow((tempXm-tempXn),2) + (float)pow((tempYm-tempYn),2) ));
	return temp_V;
	temp_V=0; tempXm=0; tempYm=0; tempXn=0; tempYn=0;
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