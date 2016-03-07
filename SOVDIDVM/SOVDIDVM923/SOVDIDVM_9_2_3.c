#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.1415926535
#define e 2.718
#define PN 1000
#define iter 10
#define GAM 1.0

float cacu_ri(float, float);
float cacu_cetai(float);

float cacu_randomNum(void);
float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);						//trans radian to degree
float cacu_9_2_1_randomNum(float);

main()
{
	float xiN[PN]={0,}, yiN[PN]={0,}, xi[PN]={0,}, yi[PN]={0,};
	float del_cetai[PN]={0,}, del_ri[PN]={0,}, del_t, Pi[PN]={0,}, randomQ[PN]={0,}, ri[PN]={0,};
	int i, t;

	FILE *Dif = fopen("Dif.dat","wt");
	del_t = 0.1;
	randomQ[0] = 0.4;
	srand((unsigned)time(NULL));

	for(t=0;t<iter;t++)
	{
		for(i=0;i<PN;i++)
		{
		//	randomQ[i] = cacu_9_2_1_randomNum(randomQ[i]);
		//	randomQ[i+1] = randomQ[i];
			randomQ[i] = cacu_randomNum();
			ri[i] = cacu_randomNum();
			printf("%0.4f\n", randomQ[i]);
		}

		for(i=0;i<PN;i++)
		{
			del_cetai[i] = cacu_cetai(randomQ[i]);
			del_ri[i] = cacu_ri(ri[i], del_t);

			xiN[i] = xi[i] + del_ri[i] * (float)cos(del_cetai[i]);
			yiN[i] = yi[i] + del_ri[i] * (float)sin(del_cetai[i]);
	//		fprintf(Dif, "%0.4f %0.4f\n", xiN[i], yiN[i]);

			xi[i] = xiN[i];
			yi[i] = yiN[i];
		}
		if(t==iter-1)
		{
			for(i=0;i<PN;i++)
			{
				fprintf(Dif, "%0.4f %0.4f\n", xiN[i], yiN[i]);
			}
		}
	}

	return 0;
}

float cacu_ri(float R, float T)
{
	float temp_ri, temp_r, v, temp_t, Pi;
	temp_r=R; temp_t=T; v=1.0;
	Pi = 1 - (float)pow(e, (((-1)*temp_r*temp_r)/(4*v*temp_t)));
	if(Pi==0)
	{
		return 0;
	}
	else
	{
		temp_ri = (float)sqrt(4*v*temp_t*(float)log(1/Pi));
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

float cacu_randomNum(void)
{
	float temp_q;
//	srand((unsigned)time(NULL));
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

float cacu_9_2_1_randomNum(float P)
{
	float temp_p, C;
	temp_p=P;
//	srand((unsigned)time(NULL));
	C=(float)(rand()%100)/100;
//	temp_p = (float)pow((1.01316+temp_p),5);
	temp_p = (float)pow((C+temp_p),5);
	temp_p = temp_p - (int)temp_p;
	return temp_p;
}
