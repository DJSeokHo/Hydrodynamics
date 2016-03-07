#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.1415926535
#define e 2.718281828459
#define GAM 1.0
#define M 50
#define GeoLength 1.0
#define Uinf 1.0
#define iter 250
#define Con M * iter

float init_cp(float, float);                    //initialization  Control points of the circle

float cacu_gamsi(float);
float cacu_gamsi_sig(float, float, float, float);

float cacu_ri(float, float);
float cacu_cetai(float);
float cacu_randomNum(void);

float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);						//trans radian to degree
float cacu_delta_uj_for_9_4_3(float, float, float, float, float);				//for convection
float cacu_delta_vj_for_9_4_3(float, float, float, float, float);				//for convection
float cacu_delta_ui_for_9_4_3_mirror(float, float);								//for mirror convection
float cacu_delta_vi_for_9_4_3_mirror(void);										//for mirror convection

main()
{
	FILE *geo=fopen("geo.dat","wt");
	FILE *cp=fopen("cp.dat","wt");
	FILE *Cov=fopen("Cov.dat","wt");
	FILE *Gam=fopen("Gam.dat","wt");

	int i, j, t, k, d, l, n;
	float del_S[M], gam_si[M], ui,  xi, xj, yj, delta_gamj[M * iter], delta_gami[M]; 
	float SIGMA[M]={0,}, Ui[M]={0,};
	float geo_x[M+1], geo_y[M+1], xmid[M], ymid[M];
	float vortices_X[M]={0,}, vortices_Y[M]={0,}, randomQ[M], ri[M], del_cetai[M], del_ri[M], deltaT;
	float randomQ_s[Con], ri_s[Con], del_cetai_s[Con], del_ri_s[Con];
	float xiN[M], yiN[M];
	float Uji[M]={0,}, Vji[M]={0,}, u_j[Con]={0,}, v_j[Con]={0,}, u_i[M]={0,}, v_i[M]={0,};
	float created_X_new[M], created_Y_new[M];
	float created_X_old[Con], created_Y_old[Con];

	deltaT = 0.015;
	srand((unsigned)time(NULL));
//************************************** set body ******************************************

	for(i=0; i<=M; i++)
	{
		geo_x[i] = (float)i/M;
		geo_y[i] = 0.0;
		fprintf(geo, "%0.4f %0.4f\n", geo_x[i], geo_y[i]);
	}
	
	for(i=0; i<M; i++)
	{
		xmid[i] = init_cp(geo_x[i+1], geo_x[i]);
		ymid[i] = init_cp(geo_y[i+1], geo_y[i]);
		del_S[i] = GeoLength / M;							//get delta_S
		fprintf(cp, "%0.4f %0.4f\n", xmid[i], ymid[i]);	
		vortices_X[i] = xmid[i];
		vortices_Y[i] = ymid[i];
	}
//******************************************************************************************
//------------------------------------------------------------------------------------------

//********************* Start Vortex creat diffusion and convection ************************

	for(t=0; t<iter; t++)
	{	

//************************************* creat vortex ***************************************
		if(t==0)
		{
			for(i=0; i<M; i++)
			{
				gam_si[i] = Uinf;									//get GAM_si				
				delta_gami[i] = gam_si[i] * del_S[i];				//new discrete vortices
				delta_gamj[i] = delta_gami[i];
			//	printf("%0.4f\n", delta_gami[i]);
				fprintf(Gam, "%0.4f\n", delta_gamj[i]);
			}
		}

		if(t>0)
		{
			for(i=0; i<M; i++)
			{		
				for(j=0; j<t*M; j++)
				{
					SIGMA[i] += cacu_gamsi_sig(vortices_X[i], created_X_old[j], created_Y_old[j], delta_gamj[j]);	//(xi, xj, yj, delta_GAM_j)
				}
				Ui[i] = cacu_gamsi(SIGMA[i]);
				delta_gami[i] = Ui[i] * del_S[i];					//new discrete vortices
			//	delta_gami[i] = Uinf * del_S[i];
				delta_gamj[i+t*M] = delta_gami[i];
				
				fprintf(Gam, "%0.4f\n", delta_gamj[i+t*M]);
			//	printf("%0.4f\n", delta_gamj[i+t*M]);
			}
			for(i=0; i<M; i++)
			{
				SIGMA[i] = 0.0;
			}
		}

//******************************************************************************************
//------------------------------------------------------------------------------------------
		
//****************************************** diffusion *************************************

		if(t>0)
		{
			for(i=0;i<t*M;i++)
			{
				randomQ_s[i] = cacu_randomNum();
				ri_s[i] = cacu_randomNum();
			//	printf("%0.4f %0.4f\n", randomQ_s[i], ri_s[i]);
			}

			for(i=0; i<t*M; i++)
			{
				del_cetai_s[i] = cacu_cetai(randomQ_s[i]);
				del_ri_s[i] = cacu_ri(ri_s[i], deltaT);
			//	printf("%0.4f %0.4f\n", del_cetai_s[i], del_ri_s[i]);

				created_X_old[i] = created_X_old[i] + del_ri_s[i] * (float)cos(del_cetai_s[i]);
				created_Y_old[i] = (float)fabs(created_Y_old[i] + del_ri_s[i] * (float)sin(del_cetai_s[i]));		

			//	printf("%0.4f %0.4f\n", created_X_old[i], created_Y_old[i]);
			}
//******************************************************************************************
//------------------------------------------------------------------------------------------


//***************************************** convection *************************************
			for(j=0;j<M*t;j++)
			{
				for(i=0;i<M;i++)
				{
					Uji[i] = cacu_delta_uj_for_9_4_3(created_X_old[j], created_Y_old[j], created_X_new[i], created_X_new[i], delta_gami[i]);		//xj, yj, xi, yi, gam
					Vji[i] = cacu_delta_vj_for_9_4_3(created_X_old[j], created_Y_old[j], created_X_new[i], created_X_new[i], delta_gami[i]);		//xj, yj, xi, yi, gam
					u_j[j] += Uji[i];
					v_j[j] += Vji[i];	
				}
			//	printf("%0.4f %0.4f\n", u_j[j], v_j[j]);
			//	printf("%0.4f %0.4f\n", (u_j[j]+Uinf) * deltaT, v_j[j]);

				created_X_old[j] = created_X_old[j] + (u_j[j] + Uinf) * deltaT;
				created_Y_old[j] = (float)fabs(created_Y_old[j] + v_j[j] * deltaT);
				
			//	printf("%0.4f %0.4f\n", created_X_old[j], created_Y_old[j]);
			}
			for(j=0;j<M*t;j++)
			{
				if(t==iter-1)
				{
					fprintf(Cov, "%0.4f %0.4f\n", created_X_old[j], created_Y_old[j]);
				}
				u_j[j] = 0.0;
				v_j[j] = 0.0;
			}
			
//******************************************************************************************
//------------------------------------------------------------------------------------------
		}
		

//****************************************** diffusion *************************************
		for(i=0;i<M;i++)
		{
			randomQ[i] = cacu_randomNum();
			ri[i] = cacu_randomNum();
		//	printf("%0.4f  %0.4f\n", randomQ[i], ri[i]);
		}

		for(i=0; i<M; i++)
		{
			del_cetai[i] = cacu_cetai(randomQ[i]);
			del_ri[i] = cacu_ri(ri[i], deltaT);

		//	printf("%0.4f  %0.4f\n", del_cetai[i], del_ri[i]);

			xiN[i] = vortices_X[i] + del_ri[i] * (float)cos(del_cetai[i]);
			yiN[i] = (float)fabs(vortices_Y[i] + del_ri[i] * (float)sin(del_cetai[i]));

		//	printf("%0.4f  %0.4f\n", xiN[i], yiN[i]);
		}
//******************************************************************************************
//------------------------------------------------------------------------------------------
		

//***************************************** convection *************************************

		for(i=0; i<M; i++)
		{
			u_i[i] = cacu_delta_ui_for_9_4_3_mirror(yiN[i], delta_gami[i]);
			v_i[i] = cacu_delta_vi_for_9_4_3_mirror();
		//	printf("%0.4f  %0.4f\n", u_i[i], v_i[i]);
		}

		for(i=0; i<M; i++)
		{
			xiN[i] = xiN[i] + (u_i[i] + Uinf) * deltaT;
			yiN[i] = (float)fabs(yiN[i] + v_i[i] * deltaT);

		//	printf("%0.4f  %0.4f\n", xiN[i], yiN[i]);

			created_X_new[i] = xiN[i];
			created_Y_new[i] = yiN[i];
		}

		for(i=0; i<M; i++)
		{
			if(t==iter-1)
			{
				fprintf(Cov, "%0.4f %0.4f\n", created_X_new[i], created_Y_new[i]);
			}

			created_X_old[i+t*M] = created_X_new[i];
			created_Y_old[i+t*M] = created_Y_new[i];
		}
		
//******************************************************************************************
//------------------------------------------------------------------------------------------


//**************************************** combination *************************************



//******************************************************************************************
//------------------------------------------------------------------------------------------
	}


	return 0;
}

float init_cp(float V_1, float V)
{
	float temp_cp_, temp_v1_, temp_v_;
	temp_v1_ = V_1; temp_v_ = V;
	temp_cp_=0.5*(temp_v1_ + temp_v_);
	return temp_cp_;
}

float cacu_gamsi(float SIGMA)
{
	float temp_gamsi, temp_sigma;
	temp_sigma=SIGMA;
	temp_gamsi = Uinf - (1 / PI) * temp_sigma;
	return temp_gamsi;
}

float cacu_gamsi_sig(float XI, float XJ, float YJ, float DELGAMJ)
{
	float temp_sig, temp_xi, temp_xj, temp_yj, temp_delgamj;
	temp_xi=XI; temp_xj=XJ; temp_yj=YJ; temp_delgamj=DELGAMJ;
	if(((float)pow((temp_xi - temp_xj),2) + (float)pow(temp_yj,2))==0.0)
	{
		temp_sig = 0.0;
	}
	else if(temp_yj==0.0)
	{
		temp_sig = 0.0;
	}
	else
	{
		temp_sig = (temp_delgamj * temp_yj) / ((float)pow((temp_xi - temp_xj),2) + (float)pow(temp_yj,2));
	}
	return temp_sig;	
}

float cacu_randomNum(void)
{
	float temp_q;
	temp_q = (float)(rand()%100)/100;
	return temp_q;
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

float cacu_ri(float R, float T)
{
	float temp_ri, temp_r, temp_v, temp_t, Pi;
	temp_r=R; temp_t=T; temp_v=1.0;
	Pi = 1 - (float)pow(e, (((-1)*temp_r*temp_r)/(4*temp_v*temp_t)));
	if(Pi==0)
	{
		return 0;
	}
	else
	{
		temp_ri = (float)sqrt(4*temp_v*temp_t*(float)log(1/Pi));
		return temp_ri;
	}
}

float cacu_cetai(float Q)
{
	float temp_ceta, temp_q;
	temp_q=Q;
	temp_ceta = 2 * PI * temp_q;
	return temp_ceta;
}

float cacu_delta_uj_for_9_4_3(float XJ, float YJ, float XI, float YI, float DELTAGAMI)
{
	float temp_uj, temp_xj, temp_yj, temp_xi, temp_yi, temp_del_gami, temp_r1_2, temp_r2_2;
	temp_xj=XJ; temp_yj=YJ; temp_xi=XI; temp_yi=YI; temp_del_gami=DELTAGAMI;
	temp_r1_2 = (float)pow( (temp_xj - temp_xi), 2) + (float)pow( (temp_yj - temp_yi), 2);
	temp_r2_2 = (float)pow( (temp_xj - temp_xi), 2) + (float)pow( (temp_yj + temp_yi), 2);
	if(temp_r1_2==0.0)
	{
	//	temp_r1_2 = 0.000001;
		temp_uj = (temp_del_gami/(2*PI)) * (((temp_yj-temp_yi))-((temp_yj+temp_yi)/temp_r2_2));
	}
	if(temp_r2_2==0.0)
	{
	//	temp_r2_2 = 0.000001;
		temp_uj = (temp_del_gami/(2*PI)) * (((temp_yj-temp_yi)/temp_r1_2)-((temp_yj+temp_yi)));
	}
	temp_uj = (temp_del_gami/(2*PI)) * (((temp_yj-temp_yi)/temp_r1_2)-((temp_yj+temp_yi)/temp_r2_2));
	return temp_uj;
}

float cacu_delta_vj_for_9_4_3(float XJ, float YJ, float XI, float YI, float DELTAGAMI)
{
	float temp_vj, temp_xj, temp_yj, temp_xi, temp_yi, temp_del_gami, temp_r1_2, temp_r2_2;
	temp_xj=XJ; temp_yj=YJ; temp_xi=XI; temp_yi=YI; temp_del_gami=DELTAGAMI;
	temp_r1_2 = (float)pow( (temp_xj - temp_xi), 2) + (float)pow( (temp_yj - temp_yi), 2);
	temp_r2_2 = (float)pow( (temp_xj - temp_xi), 2) + (float)pow( (temp_yj + temp_yi), 2);
	if(temp_r1_2==0.0)
	{
	//	temp_r1_2 = 0.00001;
		temp_vj = (temp_del_gami/(2*PI)) * (temp_xj-temp_xi) * ((1/temp_r2_2)-(1));
	}
	if(temp_r2_2==0.0)
	{
	//	temp_r2_2 = 0.00001;
		temp_vj = (temp_del_gami/(2*PI)) * (temp_xj-temp_xi) * ((1)-(1/temp_r1_2));
	}
	temp_vj = (temp_del_gami/(2*PI)) * (temp_xj-temp_xi) * ((1/temp_r2_2)-(1/temp_r1_2));
	return temp_vj;
}

float cacu_delta_ui_for_9_4_3_mirror(float YI, float DELTAGAMI)
{
	float temp_ui, temp_yi, temp_del_gami;
	temp_yi = YI; temp_del_gami=DELTAGAMI;
	if(temp_yi<0.01)
	{
		temp_ui = ((-1) * temp_del_gami);
	}
	else
	{
		temp_ui = ((-1) * temp_del_gami)/(4*PI*temp_yi);
	}
	return temp_ui;
}

float cacu_delta_vi_for_9_4_3_mirror(void)
{
	return 0;
}
