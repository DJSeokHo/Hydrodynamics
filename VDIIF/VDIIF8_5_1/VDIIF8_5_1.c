#include <stdio.h>
#include <math.h>

#define N 200
#define PI 3.1415926535

float cacu_deltaGAM(float gam, float deltaT);
float cacu_epsl(float gam, float deltaT);
float rad_to_ang(float);
float ang_to_rad(float);

main()
{
	FILE *wedge = fopen("wedge.dat", "wt");
	FILE *gamA = fopen("gamA.dat", "wt");
	FILE *gamB = fopen("gamB.dat", "wt");

	float wedgeX[N]={0,}, wedgeY[N]={0,};
	float deltaGamA, deltaGamB;
	float epslA, epslB;
	float deltaT;
	float ceta;
	float point_A_x, point_A_y, point_B_x, point_B_y;
	float epsl_length;
	int i;

	wedgeX[0]=0.0, wedgeY[0]=0.0;
	wedgeX[1]=0.6, wedgeY[1]=0.36;
	wedgeX[2]=0.6, wedgeY[2]=-0.36;
	wedgeX[3]=0.0, wedgeY[3]=0.0;
	
	epsl_length=0.025;
	ceta=30;

	point_A_x=epsl_length*(float)cos(ang_to_rad(ceta))+wedgeX[1];
	point_A_y=epsl_length*(float)sin(ang_to_rad(ceta))+wedgeY[1];

	fprintf(gamA, "%0.4f %0.4f\n", point_A_x, point_A_y);

	point_B_x=epsl_length*(float)cos(ang_to_rad(-1*ceta))+wedgeX[2];
	point_B_y=epsl_length*(float)sin(ang_to_rad(-1*ceta))+wedgeY[2];

	fprintf(gamB, "%0.4f %0.4f\n", point_B_x, point_B_y);

	for(i=0; i<4; i++)
	{
		fprintf(wedge, "%0.4f %0.4f\n", wedgeX[i], wedgeY[i]);
	}
	
	
	return 0;
}

float cacu_deltaGAM(float gam, float deltaT)
{
	float temp_deltaGAM_, temp_gam_, temp_deltaT_;
	temp_gam_=gam; temp_deltaT_=deltaT;
	temp_deltaGAM_=0.5 * (float)pow(temp_gam_,2)*temp_deltaT_;
	return temp_deltaGAM_;
	temp_deltaGAM_=0; temp_gam_=0; temp_deltaT_=0;
}

float cacu_epsl(float gam, float deltaT)
{
	float temp_epsl_, temp_gam_, temp_deltaT_;
	temp_gam_=gam; temp_deltaT_=deltaT;
	temp_epsl_=0.25 * temp_gam_*temp_deltaT_;
	return temp_epsl_;
	temp_epsl_=0; temp_gam_=0; temp_deltaT_=0;
}

float rad_to_ang(float RAD)
{
	float temp_ang_, temp_rad_;
	temp_rad_=RAD;
	temp_ang_ = temp_rad_/PI*180;
	return temp_ang_;
	temp_rad_=0; temp_ang_=0;
}

float ang_to_rad(float ANG)
{
	float temp_rad, temp_ang;
	temp_ang=ANG;
	temp_rad=(temp_ang*PI)/180;
	return temp_rad;
	temp_rad=0; temp_ang=0;
}