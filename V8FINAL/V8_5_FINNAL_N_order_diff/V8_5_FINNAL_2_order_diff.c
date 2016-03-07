#include <stdio.h>
#include <math.h>

#define M 11					//number of surfaces
#define N 2000				//Array size
#define PI 3.1415926535		//PI
#define a 0.5				//R of cycle
#define Winf 1.0			//Winf
#define flow_angle 0.0		//flow angle
#define iter 200				//number of iterations
#define difference 2

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

float cacu_Umn(float, float, float, float);		//Get Umn. Here m & n means convection shed points only 
float cacu_Vmn(float, float, float, float);		//Get Vmn. Here m & n means convection shed points only 

main()
{
	draw();
	return 0;
}

void draw(void)
{
	int i, j, k, time, VortexNum, VortexIndex, d;		//for Loop
	FILE *circle=fopen("circle.dat","wt");				//for geometry
	FILE *pivotalP=fopen("pivotalP.dat","wt");			//for geometry
	FILE *driftPathA=fopen("driftPathA.dat","wt");		//for driftPath A
	FILE *driftPathB=fopen("driftPathB.dat","wt");		//for driftPath B
	FILE *shedA=fopen("shedA.dat","wt");				//for shed points A
	FILE *shedB=fopen("shedB.dat","wt");				//for shed points B

	float cycx[N]={0,}, cycy[N]={0,}, xmid[N]={0,}, ymid[N]={0,}, deltaSn[N]={0,}, angle, slope[N]={0,}, cosbetam[N]={0,}, sinbetam[N]={0,};		//init the cycle
	
	float Uinf, Vinf;	
	float Coeff[N]={0,};																	//for Martensen
	float K_Sm_Sn_Winf[N][N]={0,}, inv_Kmn_Winf[N][N]={0,}, Surface_gam_Winf[N]={0,}, rhsM_Winf[N]={0,};	//for Martensen

	float K_Sm_SnA[N][N]={0,}, inv_KmnA[N][N]={0,}, Surface_gamA[N]={0,}, rhsMA[N]={0,};	//for Martensen
	float Surface_gam_addA[N]={0,};															//for Martensen
	
	float K_Sm_SnB[N][N]={0,}, inv_KmnB[N][N]={0,}, Surface_gamB[N]={0,}, rhsMB[N]={0,};	//for Martensen
	float Surface_gam_addB[N]={0,};															//for Martensen

	float Surface_gamAB[N]={0,};

	float rmjA[N]={0,}, UmjA[N]={0,}, VmjA[N]={0,};
	float UmnA[N]={0,}, VmnA[N]={0,};

	float rmjB[N]={0,}, UmjB[N]={0,}, VmjB[N]={0,};
	float UmnB[N]={0,}, VmnB[N]={0,};

	float delta_GAM_A, delta_GAM_B, deltaT, epsl_length, ceta;							//two shed points
	float delta_GAM_AB[N];

	float shed_A_x, shed_A_y;															//two shed points
	float shed_B_x, shed_B_y;															//two shed points
	float point_A_x[N]={0,}, point_A_y[N]={0,};											//two shed points
	float point_B_x[N]={0,}, point_B_y[N]={0,};											//two shed points
	float point_AB_x[N]={0,}, point_AB_y[N]={0,};										//two shed points
	

	float point_A_x_next[N]={0,}, point_A_y_next[N]={0,};								//for drift path
	float point_B_x_next[N]={0,}, point_B_y_next[N]={0,};								//for drift path

	float sigma_rhsm_for_mjA[N][N]={0,}, sigma_rhsm_mjA[N]={0,};						//for K function's second sigma item
	float sigma_rhsm_for_mjB[N][N]={0,}, sigma_rhsm_mjB[N]={0,};

	float udmA[N]={0,}, vdmA[N]={0,};
	float udm_sig_oneA[N]={0,}, udm_sig_twoA[N]={0,}, vdm_sig_oneA[N]={0,}, vdm_sig_twoA[N]={0,};

	float udmB[N]={0,}, vdmB[N]={0,};
	float udm_sig_oneB[N]={0,}, udm_sig_twoB[N]={0,}, vdm_sig_oneB[N]={0,}, vdm_sig_twoB[N]={0,};

	float udmAB[N]={0,}, vdmAB[N]={0,};
	float udm_sig_oneAB[N]={0,}, udm_sig_twoAB[N]={0,}, vdm_sig_oneAB[N]={0,}, vdm_sig_twoAB[N]={0,};

	float outputstreamL_A_x[N]={0,}, outputstreamL_A_y[N]={0,};
	float outputstreamL_B_x[N]={0,}, outputstreamL_B_y[N]={0,};

	float temp_u_A[N], temp_v_A[N], temp_x_A[N], temp_y_A[N], temp_ever_u_A[N][N], temp_ever_v_A[N][N];		//for N-order difference of A
	float temp_u_B[N], temp_v_B[N], temp_x_B[N], temp_y_B[N], temp_ever_u_B[N][N], temp_ever_v_B[N][N];		//for N-order difference of B
	float temp_sam_u_A[N]={0,}, temp_sam_v_A[N]={0,}, temp_sam_u_B[N]={0,}, temp_sam_v_B[N]={0,};

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
		}
		rhsM_Winf[i] = cacu_rhsm_Winf(Uinf, cosbetam[i], Vinf, sinbetam[i]);	//(U, cos, V, sin)
	}
	matrix_inversion_using_elementary_operation_Winf(K_Sm_Sn_Winf,inv_Kmn_Winf,M-1);
	for(i=0;i<M-1;i++)
	{
		Surface_gam_Winf[i]=0.0;
		for(j=0;j<M-1;j++)
		{
			Surface_gam_Winf[i]+=inv_Kmn_Winf[i][j]*rhsM_Winf[j];				//Get surface GAM
		}
	}

	deltaT = 0.05;																	//T is using for make shed point
	
//*************************************||            Shed two free vortex           ||*************************************

//	epsl_length = cacu_epsl(deltaSn[2]);											//init shed points
//	ceta = 45;	

//	shed_A_x=epsl_length*(float)cos(ang_to_rad(ceta))+xmid[2];
//	shed_A_y=epsl_length*(float)sin(ang_to_rad(ceta))+ymid[2];

//	shed_B_x=epsl_length*(float)cos(ang_to_rad(ceta))+xmid[7];
//	shed_B_y=ymid[7]-epsl_length*(float)sin(ang_to_rad(ceta));

	shed_A_x=xmid[1];
	shed_A_y=ymid[1];

	shed_B_x=xmid[2];
	shed_B_y=ymid[2];

	fprintf(shedA, "%0.4f %0.4f\n", shed_A_x, shed_A_y);
	fprintf(shedB, "%0.4f %0.4f\n", shed_B_x, shed_B_y);

	outputstreamL_A_x[iter+1] = shed_A_x;
	outputstreamL_A_y[iter+1] = shed_A_y;
	outputstreamL_B_x[iter+1] = shed_B_x;
	outputstreamL_B_y[iter+1] = shed_B_y;

	delta_GAM_A = cacu_deltaGAM(Surface_gam_Winf[0], deltaT);				//  0.0538
	delta_GAM_B = cacu_deltaGAM(Surface_gam_Winf[1], deltaT);				// -0.0538

//	delta_GAM_A = 0.25;
//	delta_GAM_B = 0.25;

	for(time = 0; time < iter; time++)		//time step
	{
		VortexNum = time + 1;			//Vortex number follow time steps

		point_A_x[time] = shed_A_x;
		point_A_y[time] = shed_A_y;

		point_B_x[time] = shed_B_x;
		point_B_y[time] = shed_B_y;




		for(d=0;d<difference;d++)		//for N-order difference
		{








//*************************************||           Caculate the Umj & Vmj          ||*************************************
	
			for(i=0;i<VortexNum;i++)				//Vortices index
			{
				for(j=0;j<M-1;j++)	//Here means only have one free vortex
				{
					rmjA[j] = cacu_rm(xmid[j], point_A_x[i], ymid[j], point_A_y[i]);
					UmjA[j] = cacu_Um(ymid[j], point_A_y[i], rmjA[j]);
					VmjA[j] = cacu_Vm(xmid[j], point_A_x[i], rmjA[j]);

					sigma_rhsm_for_mjA[i][j] = cacu_rhsm(UmjA[j], cosbetam[j], VmjA[j], sinbetam[j], delta_GAM_A);

					rmjB[j] = cacu_rm(xmid[j], point_B_x[i], ymid[j], point_B_y[i]);
					UmjB[j] = cacu_Um(ymid[j], point_B_y[i], rmjB[j]);
					VmjB[j] = cacu_Vm(xmid[j], point_B_x[i], rmjB[j]);
					
					sigma_rhsm_for_mjB[i][j] = cacu_rhsm(UmjB[j], cosbetam[j], VmjB[j], sinbetam[j], delta_GAM_B);	
				}
			}

//*************************************|| Martensen analysis Umj & Vmj of shed point||*************************************

			for(i=0;i<M-1;i++)
			{
				sigma_rhsm_mjA[i] = 0.0;
				sigma_rhsm_mjB[i] = 0.0;
				for(j=0;j<VortexNum;j++)
				{
					sigma_rhsm_mjA[i] += sigma_rhsm_for_mjA[j][i];
					sigma_rhsm_mjB[i] += sigma_rhsm_for_mjB[j][i];
				}
			}		
			for(i=0;i<M-1;i++)
			{
				rhsMA[i] = rhsM_Winf[i] - sigma_rhsm_mjA[i];
				rhsMB[i] = rhsM_Winf[i] - sigma_rhsm_mjB[i];
			}			
			for(i=0;i<M-1;i++)																//Martensen analysis
			{
				for(j=0;j<M-1;j++)
				{
					K_Sm_SnA[i][j] = cacu_K_Sm_Sn(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);	//(XM, YM, XN, YN, COS, SIN)
					K_Sm_SnB[i][j] = cacu_K_Sm_Sn(Coeff[j], xmid[i],ymid[i],xmid[j],ymid[j],cosbetam[i],sinbetam[i]);	//(XM, YM, XN, YN, COS, SIN)
				}
			}	
			matrix_inversion_using_elementary_operation(K_Sm_SnA,inv_KmnA,M-1);
			for(i=0;i<M-1;i++)
			{
				Surface_gamA[i]=0.0;
				for(j=0;j<M-1;j++)
				{
					Surface_gamA[i]+=inv_KmnA[i][j]*rhsMA[j];				//Get surface GAM
				}
			}
			matrix_inversion_using_elementary_operation(K_Sm_SnB,inv_KmnB,M-1);
			for(i=0;i<M-1;i++)
			{
				Surface_gamB[i]=0.0;
				for(j=0;j<M-1;j++)
				{
					Surface_gamB[i]+=inv_KmnB[i][j]*rhsMB[j];				//Get surface GAM
				}
			}
			for(i=0;i<M-1;i++)						//Plus all surface gam into a array
			{
				Surface_gamAB[i] = Surface_gamA[i] + Surface_gamB[i];
			}

//*************************************||     Cacu the convection of Umn & Vmn      ||*************************************

			for(i=0;i<VortexNum;i++)				//Plus all points into a array
			{
				point_AB_x[i] = point_A_x[i];
				point_AB_y[i] = point_A_y[i];
			}
			for(i=0;i<VortexNum;i++)				//Plus all points into a array
			{
				point_AB_x[i+VortexNum] = point_B_x[i];
				point_AB_y[i+VortexNum] = point_B_y[i];
			}

			for(i=0;i<VortexNum;i++)					//Plus gam A&B into a array
			{
				delta_GAM_AB[i] = delta_GAM_A;
			}
			for(i=0;i<VortexNum;i++)					//Plus gam A&B into a array
			{
				delta_GAM_AB[i+VortexNum] = delta_GAM_B;
			}

			for(i=0;i<VortexNum*2;i++)
			{
				udm_sig_oneAB[i] = 0.0;
				vdm_sig_oneAB[i] = 0.0;
				for(j=0;j<VortexNum*2;j++)
				{
					if(i!=j)
					{
						udm_sig_oneAB[i] += delta_GAM_AB[i] * cacu_Umn(point_AB_x[i], point_AB_y[i], point_AB_x[j], point_AB_y[j]);
						vdm_sig_oneAB[i] += delta_GAM_AB[i] * cacu_Vmn(point_AB_x[i], point_AB_y[i], point_AB_x[j], point_AB_y[j]);
					}
				}
			}
		
			for(i=0;i<VortexNum;i++)									//Distribution sig_one into two arrays
			{
				udm_sig_oneA[i] = udm_sig_oneAB[i];
				vdm_sig_oneA[i] = vdm_sig_oneAB[i];
			}
			for(i=0;i<VortexNum;i++)									//Distribution sig_one into two arrays
			{
				udm_sig_oneB[i] = udm_sig_oneAB[i+VortexNum];
				vdm_sig_oneB[i] = vdm_sig_oneAB[i+VortexNum];
			}
			
			for(i=0;i<VortexNum;i++)									
			{
				udm_sig_twoA[i] = 0.0; vdm_sig_twoA[i] = 0.0;
				udm_sig_twoB[i] = 0.0; vdm_sig_twoB[i] = 0.0;
				for(j=0;j<M-1;j++)
				{
					udm_sig_twoA[i] += Surface_gamA[j] * deltaSn[j] * UmjA[j];
					vdm_sig_twoA[i] += Surface_gamA[j] * deltaSn[j] * VmjA[j];
					udm_sig_twoB[i] += Surface_gamB[j] * deltaSn[j] * UmjB[j];
					vdm_sig_twoB[i] += Surface_gamB[j] * deltaSn[j] * VmjB[j];
				}
			}

//*************************************||   prepare velocity for Vortex convection  ||*************************************
			
			for(i=0;i<VortexNum;i++)									//Add two item
			{
				if(time<10)
				{
					Vinf = 0.50;
				}
				else
				{
					Vinf = 0.0;
				}

				udmA[i] = udm_sig_oneA[i] + udm_sig_twoA[i] + Uinf;
				vdmA[i] = vdm_sig_oneA[i] + vdm_sig_twoA[i] + Vinf;
				udmB[i] = udm_sig_oneB[i] + udm_sig_twoB[i] + Uinf;
				vdmB[i] = vdm_sig_oneB[i] + vdm_sig_twoB[i] + Vinf;
			}

			for(i=0;i<VortexNum;i++)									//save temp u&v
			{
				temp_u_A[i] = udmA[i];
				temp_v_A[i] = vdmA[i];
				temp_u_B[i] = udmB[i];
				temp_v_B[i] = vdmB[i];
			}
			
			for(i=0;i<VortexNum;i++)									//save temp u&v forever
			{
				temp_ever_u_A[d][i] = udmA[i];
				temp_ever_v_A[d][i] = vdmA[i];
				temp_ever_u_B[d][i] = udmB[i];
				temp_ever_v_B[d][i] = vdmB[i];
			}

//*************************************||              Vortex displace              ||*************************************
			
			for(i=0;i<VortexNum;i++)									//save displacement of A&B		
			{
				temp_x_A[i] = point_A_x[i] + temp_u_A[i] * deltaT;
				temp_y_A[i] = point_A_y[i] + temp_v_A[i] * deltaT;
				temp_x_B[i] = point_B_x[i] + temp_u_B[i] * deltaT;
				temp_y_B[i] = point_B_y[i] + temp_v_B[i] * deltaT;
			}

			for(i=0;i<VortexNum;i++)									//copy displacement for new A&B		
			{
				point_A_x[i] = temp_x_A[i];
				point_A_y[i] = temp_y_A[i];
				point_B_x[i] = temp_x_B[i];
				point_B_y[i] = temp_y_B[i];
			}

			
			

		}		//for N-order difference
		
		for(i=0;i<VortexNum;i++)									
		{
			temp_sam_u_A[i] = 0.0;
			temp_sam_v_A[i] = 0.0;
			temp_sam_u_B[i] = 0.0;
			temp_sam_v_B[i] = 0.0;
			for(d=0;d<difference;d++)
			{
				temp_sam_u_A[i] += temp_ever_u_A[d][i];
				temp_sam_v_A[i] += temp_ever_v_A[d][i];
				temp_sam_u_B[i] += temp_ever_u_B[d][i];
				temp_sam_v_B[i] += temp_ever_v_B[d][i];
			}
		}

		for(i=0;i<VortexNum;i++)									
		{
			point_A_x_next[i] = point_A_x[i] + (1/difference) * temp_sam_u_A[i] * deltaT;
			point_A_y_next[i] = point_A_y[i] + (1/difference) * temp_sam_v_A[i] * deltaT;
			point_B_x_next[i] = point_B_x[i] + (1/difference) * temp_sam_u_B[i] * deltaT;
			point_B_y_next[i] = point_B_y[i] + (1/difference) * temp_sam_v_B[i] * deltaT;
		}



//*************************************||             Output streamline             ||*************************************

		if(time == iter - 1)
		{
			for(i=0;i<VortexNum;i++)								
			{
				outputstreamL_A_x[i+1] = point_A_x_next[i];
				outputstreamL_A_y[i+1] = point_A_y_next[i];
				outputstreamL_B_x[i+1] = point_B_x_next[i];
				outputstreamL_B_y[i+1] = point_B_y_next[i];
			}
		}

		for(i=0;i<VortexNum;i++)									
		{
			point_A_x[i] = point_A_x_next[i];
			point_A_y[i] = point_A_y_next[i];
			point_B_x[i] = point_B_x_next[i];
			point_B_y[i] = point_B_y_next[i];
		}

		if(time == iter - 1)
		{
			for(i=VortexNum+1;i>0;i--)
			{
				fprintf(driftPathA, "%0.4f %0.4f\n", outputstreamL_A_x[i], outputstreamL_A_y[i]);
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
	temp_rhsm=temp_gam*((temp_u*temp_cos)+(temp_v*temp_sin));
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