#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.1415926535
#define e 2.718
#define PN 50

float cacu_omega(float, float);
float cacu_p(float, float, float, float);
float cacu_ri(float, float);

float cacu_randomNum(void);
float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);						//trans radian to degree
float cacu_9_2_1_randomNum(float);

main()
{
	float Omega, r, strength, p, del_ceta, del_r, randomQ, Cetai, ri;
	int i;

	FILE *ApointV = fopen("ApointV.dat","wt");
	
	del_r = 0.1;
	del_ceta = 1.0;

	for(i=-50;i<=50;i++)
	{
		r = (float)i/10;
		Omega = cacu_omega(r,1.0);
		fprintf(ApointV,"%0.4f %0.4f\n", r, Omega);
	}
	
	p = cacu_p(0.5, 0.25, del_r, del_ceta); //r,   t,   delta r,   delta ceta
//	printf("probability is  %0.4f\n", p);

	Cetai = 2 * PI * cacu_randomNum();
	ri = cacu_ri(0.1, 0.1);
	printf("Cetai is  %0.4f\n", rad_to_ang(Cetai));
	printf("Ri is  %0.4f\n", ri);

	return 0;
}

float cacu_omega(float R, float T)
{
	float GAM, v, temp_r, temp_t, temp_o;
	GAM=1.0; v=1.0; temp_r=R; temp_t=T;
	temp_o = (GAM/(4*PI*v*temp_t)) * (float)pow(e,(((-1)*temp_r*temp_r)/(4*v*temp_t)));
	return temp_o;
	temp_o=0; GAM=0; v=0; temp_r=0; temp_t=0;
}

float cacu_p(float R, float T, float delR, float delCeta)
{
	float GAM, v, temp_r, temp_t, temp_p, temp_del_r, temp_del_ceta;
	GAM=1.0/PN; v=1.0; temp_r=R; temp_t=T; temp_del_r=delR; temp_del_ceta=delCeta;
	temp_p = (1/(4*PI*v*temp_t)) * (float)pow(e,(((-1)*temp_r*temp_r)/(4*v*temp_t))) * temp_r * temp_del_r * temp_del_ceta;
	return temp_p;
	temp_p=0; GAM=0; v=0; temp_r=0; temp_t=0; temp_del_r=0; temp_del_ceta=0;
}

float cacu_ri(float R, float T)
{
	float temp_ri, temp_r, v, temp_t, Pi;
	temp_r=R; temp_t=T; v=1.0;
	Pi = 1 - (float)pow(e, (((-1)*temp_r*temp_r)/(4*v*temp_t)));
	temp_ri = (float)sqrt(4*v*temp_t*(float)log(1/Pi));
	return temp_ri;
	temp_ri=0; temp_r=0; temp_t=0; v=0;
}

float cacu_randomNum(void)
{
	float temp_q;
	srand((unsigned)time(NULL));
	temp_q = (float)(rand()%100)/100;
	return temp_q;
	temp_q=0;
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

float cacu_9_2_1_randomNum(float P)
{
	float temp_p;
	temp_p=P;
	temp_p = (float)pow((1.01316+temp_p),5);
	temp_p = temp_p - (int)temp_p;
	return temp_p;
	temp_p=0;
}
