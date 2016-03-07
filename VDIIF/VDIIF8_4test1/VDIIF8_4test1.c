#include <stdio.h>
#include <math.h>

#define PI 3.1415926535
#define N 31

#define a 0.5
#define Uinf 1.0		

float cacu_rm(float, float, float, float);      //Get Rm
float cacu_Um(float, float, float);             //Get Um
float cacu_Vm(float, float, float);             //Get Vm
float cacu_um(float, float);             //Get um
float cacu_vm(float, float);             //Get vm
float cacu_qm(float, float);                    //Get qm
float init_x(float, float);                     //initialization  X of the circle
float init_y(float, float);                     //initialization  Y of the circle
float init_cp(float, float);                    //initialization  Control points of the circle
float ang_to_rad(float);                        //trans degree to radian
float cacu_xb(float, float, float);
float cacu_yb(float, float, float);
float cacu_xd(float, float, float, float);
float cacu_yd(float, float, float, float);

main()
{
	float x[N], y[N], xm[N], ym[N], xja, yja, xjb, yjb, xjd, yjd;
	float Umj[N], Vmj[N];
	float um[N]={0,}, vm[N]={0,};
	float ub = 0, vb = 0;
	float ud = 0, vd = 0;
	float rmj[N], qm[N];
	float deltaT;
	float deltaGamSM[N];
	float deltaSm[N];
	float angle;
	float deltaGam;
	int i, time, iter;

	FILE *circle=fopen("circle.dat","wt");
	FILE *pivotalP=fopen("pivotalP.dat","wt");
	FILE *driftPath=fopen("driftPath.dat","wt");

	for(i = 0, angle = 0; i <= N-1, angle <= 360; i++, angle += (360/(N-1)))
	{
		x[i] = init_x(ang_to_rad(angle),a);
		y[i] = init_y(ang_to_rad(angle),a);
		fprintf(circle, "%0.4f  %0.4f\n", x[i], y[i] );
	}

	for(i = 0; i < N - 1; i++)
	{
		xm[i] = init_cp(x[i+1],x[i]);
		ym[i] = init_cp(y[i+1],y[i]);
		fprintf(pivotalP, "%0.4f  %0.4f\n", xm[i], ym[i] );
	}

	for(i = 0; i < N - 1; i++)
	{
		deltaSm[i] = (float)sqrt((float)pow((x[i+1] - x[i]),2) + (float)pow((y[i+1] - y[i]),2));
	}

	xja=-1.0; yja=0.2;
	deltaT = 0.05;
	deltaGam = 0.05;
	iter = 100;

	fprintf(driftPath,"%0.4f %0.4f\n", xja,yja);

//**************************************** from a to b ****************************************
	for(time=0; time<iter; time++)
	{
		for(i = 0; i < N - 1; i++)
		{
			rmj[i] = cacu_rm(xm[i], xja, ym[i], yja);
			Umj[i] = cacu_Um(ym[i], yja, rmj[i]);
			Vmj[i] = cacu_Vm(xm[i], xja, rmj[i]);

			um[i] = cacu_um(Umj[i], deltaGam);
			vm[i] = cacu_vm(Vmj[i], deltaGam);
			qm[i] = cacu_qm(um[i], vm[i]);

			deltaGamSM[i] = qm[i];

			ub += (deltaGamSM[i]*deltaSm[i]*Umj[i]);
			vb += (deltaGamSM[i]*deltaSm[i]*Vmj[i]);
		}

		ub = (-1)*ub;
		vb = (-1)*vb;

		xjb = xja + (ub) * deltaT;
		yjb = yja + (vb) * deltaT;

//**************************************** from b to d ****************************************

		for(i = 0; i < N - 1; i++)
		{

			rmj[i] = cacu_rm(xm[i], xjb, ym[i], yjb);
			Umj[i] = cacu_Um(ym[i], yjb, rmj[i]);
			Vmj[i] = cacu_Vm(xm[i], xjb, rmj[i]);
			
			um[i] = cacu_um(Umj[i], deltaGam);
			vm[i] = cacu_vm(Vmj[i], deltaGam);
			qm[i] = cacu_qm(um[i], vm[i]);

			deltaGamSM[i] = qm[i];

			ud += (deltaGamSM[i]*deltaSm[i]*Umj[i]);
			vd += (deltaGamSM[i]*deltaSm[i]*Vmj[i]);
		}

		ud = (-1)*ud;
		vd = (-1)*vd;

		xjd = xja + 0.5*(ub+ud) * deltaT;
		yjd = yja + 0.5*(vb+vd) * deltaT;
		
		fprintf(driftPath,"%0.4f %0.4f\n", xjd,yjd);
		
		xja=xjd;yja=yjd;
		ub=0;vb=0;ud=0;vd=0;
	}

	return 0;
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

float cacu_um(float Mj, float deltaGamj)
{
	float temp_um_, temp_mj_, temp_deltaGamj_;
	temp_mj_=Mj; temp_deltaGamj_=deltaGamj;
	temp_um_=temp_mj_ * temp_deltaGamj_;
	return temp_um_;
	temp_um_=0; temp_mj_=0; temp_deltaGamj_=0;
}

float cacu_vm(float Mj, float deltaGamj)
{
	float temp_vm_, temp_mj_, temp_deltaGamj_;
	temp_mj_=Mj; temp_deltaGamj_=deltaGamj;
	temp_vm_=temp_mj_ * temp_deltaGamj_;
	return temp_vm_;
	temp_vm_=0; temp_mj_=0; temp_deltaGamj_=0;
}

float cacu_rm(float XM, float X_, float YM, float Y_)
{
	float temp_rm_, temp_xm_, temp_x_, temp_ym_, temp_y_;
	temp_xm_=XM; temp_x_=X_; temp_ym_=YM; temp_y_=Y_;
	temp_rm_ = (float)sqrt( (float)pow((temp_xm_-temp_x_),2) + (float)pow((temp_ym_-temp_y_),2) );
	return temp_rm_;
	temp_xm_=0; temp_x_=0;temp_ym_=0;temp_y_=0;temp_rm_=0;
}

float cacu_qm(float UM, float VM)
{
	float temp_q_, temp_um_, temp_vm_;
	temp_um_=UM; temp_vm_=VM;
	temp_q_ = (float)sqrt( (float)pow(temp_um_,2) + (float)pow(temp_vm_,2) );
	return temp_q_;
	temp_q_=0; temp_um_=0; temp_vm_=0;
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

float cacu_xb(float Xa, float u, float T)
{
	float temp_xb, tempXa, tempu, tempT;
	tempXa=Xa; tempu=u; tempT=T;
	temp_xb=tempXa+(tempu*tempT);
	return temp_xb;
	temp_xb=0; tempXa=0; tempu=0; tempT=0;
}

float cacu_yb(float Ya, float v, float T)
{
	float temp_yb, tempYa, tempv, tempT;
	tempYa=Ya; tempv=v; tempT=T;
	temp_yb=tempYa+(tempv*tempT);
	return temp_yb;
	temp_yb=0; tempYa=0; tempv=0; tempT=0;
}

float cacu_xd(float Xa, float ua, float ub, float T)
{
	float temp_xd, tempXa, tempua, tempub, tempT;
	tempXa=Xa; tempua=ua; tempub=ub; tempT=T;
	temp_xd=tempXa+(0.5*(tempua+tempub)*tempT);
	return temp_xd;
	temp_xd=0; tempXa=0; tempua=0; tempub=0; tempT=0;
}

float cacu_yd(float Ya, float va, float vb, float T)
{
	float temp_yd, tempYa, tempva, tempvb, tempT;
	tempYa=Ya; tempva=va; tempvb=vb; tempT=T;
	temp_yd=tempYa+(0.5*(tempva+tempvb)*tempT);
	return temp_yd;
	temp_yd=0; tempYa=0; tempva=0; tempvb=0; tempT=0;
}
