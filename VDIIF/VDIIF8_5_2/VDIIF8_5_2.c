#include <stdio.h>
#include <math.h>

#define N 31
#define PI 3.1415926535
#define a 0.5


float init_x(float, float);                     
float init_y(float, float);                     
float init_cp(float, float);                    
float ang_to_rad(float);                        
float rad_to_ang(float);

float cacu_U_821(float, float, float, float);
float cacu_V_821(float, float, float, float);
float cacu_u_821(float, float);
float cacu_v_821(float, float);

float cacu_rm_842(float, float, float, float);  
float cacu_Um_842(float, float, float);             
float cacu_Vm_842(float, float, float);             
float cacu_um_842(float, float);             
float cacu_vm_842(float, float);             
float cacu_qm_842(float, float); 

void drawing(void);

main()
{
	drawing();
	return 0;
}

void drawing(void)
{
	int i, j;

	float deltaT;

	FILE *circle=fopen("circle.dat","wt");
	FILE *pivotalP=fopen("pivotalP.dat","wt");
	FILE *shedP=fopen("shedP.dat","wt");

	float angle;
	float cycx[N]={0,}, cycy[N]={0,};
	float xmid[N]={0,}, ymid[N]={0,};
	float deltaS_Element[N]={0,};

	float epsl, deltaGAM_A_angle, deltaGAM_A, deltaGAM_A_x, deltaGAM_A_y;

	float Umn_821_[N]={0,}, Vmn_821_[N]={0,};
	float u_821_=0, v_821_=0;

	float rmj_842_[N]={0,}, qm_842_[N]={0,};
	float Umj_842_[N]={0,}, Vmj_842_[N]={0,};
	float um_842_[N]={0,}, vm_842_[N]={0,};
	float u_842_=0, v_842_=0;
	float surfaceGAM_842_[N]={0,};



	//**-----init the cycle-----**
	for(i = 0, angle = 0; i <= N-1, angle <= 360; i++, angle += (360/(N-1)))
	{
		cycx[i] = init_x(ang_to_rad(angle),a);
		cycy[i] = init_y(ang_to_rad(angle),a);
		fprintf(circle, "%0.4f  %0.4f\n", cycx[i], cycy[i] );
	}



	//**-----init the middle points of per surface elements-----**
	for(i = 0; i < N - 1; i++)
	{
		xmid[i] = init_cp(cycx[i+1],cycx[i]);
		ymid[i] = init_cp(cycy[i+1],cycy[i]);
		fprintf(pivotalP, "%0.4f  %0.4f\n", xmid[i], ymid[i] );
	}
	


	//**-----init the length of each two elements-----**
	for(i = 0; i < N - 1; i++)	//element's deltaS
	{
		deltaS_Element[i] = (float)sqrt((float)pow((cycx[i+1] - cycx[i]),2) + (float)pow((cycy[i+1] - cycy[i]),2));
	//	printf("%0.4f\n", deltaS_Element[i]);
	}



	//**-----init the point delta GAM A-----**
	epsl = deltaS_Element[0]/2;
	deltaGAM_A_angle = ang_to_rad(45);
	deltaGAM_A_x = epsl * (float)cos(deltaGAM_A_angle) + xmid[7];
	deltaGAM_A_y = epsl * (float)sin(deltaGAM_A_angle) + ymid[7];

	fprintf(shedP, "%0.4f %0.4f\n", deltaGAM_A_x, deltaGAM_A_y);


	//**-----caculate the u and v of 8.2.1-----**
	for(i=0;i<1;i++)		//point delta GAM A to the 30 points in the cycle!!!
	{
		for(j=0;j<N-1;j++)
		{
		//	if(i!=j)		//because point delta GAM A is not equal any point in the cycle!!!
		//	{
				Umn_821_[j] = cacu_U_821(deltaGAM_A_x, deltaGAM_A_y, xmid[j], ymid[j]);
				Vmn_821_[j] = cacu_V_821(deltaGAM_A_x, deltaGAM_A_y, xmid[j], ymid[j]);
				u_821_ += cacu_u_821(Umn_821_[j], deltaGAM_A);
				v_821_ += cacu_v_821(Vmn_821_[j], deltaGAM_A);
		//	}
		}
	//	printf("%0.4f %0.4f\n", u_821_, v_821_);
	}
	


	//**-----caculate the u and v of 8.4.2-----**
	for(i = 0; i < N - 1; i++)
	{
		rmj_842_[i] = cacu_rm_842(xmid[i], deltaGAM_A_x, ymid[i], deltaGAM_A_y);
		Umj_842_[i] = cacu_Um_842(ymid[i], deltaGAM_A_y, rmj_842_[i]);
		Vmj_842_[i] = cacu_Vm_842(xmid[i], deltaGAM_A_x, rmj_842_[i]);
		
		um_842_[i] = cacu_um_842(Umj_842_[i], deltaGAM_A);		//point delta GAM A to the surface control points
		vm_842_[i] = cacu_vm_842(Vmj_842_[i], deltaGAM_A);		//point delta GAM A to the surface control points
		qm_842_[i] = cacu_qm_842(um_842_[i], vm_842_[i]);
		surfaceGAM_842_[i] = qm_842_[i];

		u_842_ += surfaceGAM_842_[i]*deltaS_Element[i]*Umj_842_[i];
		v_842_ += surfaceGAM_842_[i]*deltaS_Element[i]*Vmj_842_[i];
	}
//	printf("%0.4f %0.4f\n", u_821_, v_821_);
	//At the end must init u and v of 8.2.1 and 8.4.2 and the last (u,v)!!!
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

float cacu_U_821(float Xm, float Ym, float Xn, float Yn)
{
	float temp_U, tempXm, tempYm, tempXn, tempYn;
	tempXm=Xm; tempYm=Ym; tempXn=Xn; tempYn=Yn;
	temp_U=(1/(2*PI))*( (tempYm-tempYn) / ( (float)pow((tempXm-tempXn),2) + (float)pow((tempYm-tempYn),2) ));
	return temp_U;
	temp_U=0; tempXm=0; tempYm=0; tempXn=0; tempYn=0;
}

float cacu_V_821(float Xm, float Ym, float Xn, float Yn)
{
	float temp_V, tempXm, tempYm, tempXn, tempYn;
	tempXm=Xm; tempYm=Ym; tempXn=Xn; tempYn=Yn;
	temp_V=(-1/(2*PI))*((tempXm-tempXn) / ( (float)pow((tempXm-tempXn),2) + (float)pow((tempYm-tempYn),2) ));
	return temp_V;
	temp_V=0; tempXm=0; tempYm=0; tempXn=0; tempYn=0;
}
float cacu_u_821(float U, float GAM)
{
	float temp_u, tempU, tempGAM;
	tempU=U; tempGAM=GAM;
	temp_u=tempGAM * tempU;
	return temp_u;
	temp_u=0; tempU=0; tempGAM=0;
}
float cacu_v_821(float V, float GAM)
{
	float temp_v, temp, tempGAMV;
	tempV=V; tempGAM=GAM;
	temp_v=tempGAM * tempV;
	return temp_v;
	temp_v=0; tempV=0; tempGAM=0;
}

float cacu_rm_842(float XM, float X_, float YM, float Y_)
{
	float temp_rm_, temp_xm_, temp_x_, temp_ym_, temp_y_;
	temp_xm_=XM; temp_x_=X_; temp_ym_=YM; temp_y_=Y_;
	temp_rm_ = (float)sqrt( (float)pow((temp_xm_-temp_x_),2) + (float)pow((temp_ym_-temp_y_),2) );
	return temp_rm_;
	temp_xm_=0; temp_x_=0;temp_ym_=0;temp_y_=0;temp_rm_=0;
}

float cacu_Um_842(float YM, float Y_, float RM_)
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

float cacu_Vm_842(float XM, float X_, float RM_)
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

float cacu_um_842(float Mj, float Gamj)
{
	float temp_um_, temp_mj_, temp_gamj_;
	temp_mj_=Mj; temp_gamj_=Gamj;
	temp_um_=temp_mj_ * temp_gamj_;
	return temp_um_;
	temp_um_=0; temp_mj_=0; temp_gamj_=0;
}

float cacu_vm_842(float Mj, float Gamj)
{
	float temp_vm_,temp_mj_, temp_gamj_;
	temp_mj_=Mj; temp_gamj_=Gamj;
	temp_vm_=temp_mj_ * temp_gamj_;
	return temp_vm_;
	temp_vm_=0; temp_mj_=0; temp_gamj_=0;
}

float cacu_qm_842(float UM, float VM)
{
	float temp_q_, temp_um_, temp_vm_;
	temp_um_=UM; temp_vm_=VM;
	temp_q_ = (float)sqrt( (float)pow(temp_um_,2) + (float)pow(temp_vm_,2) );
	return temp_q_;
	temp_q_=0; temp_um_=0; temp_vm_=0;
}               