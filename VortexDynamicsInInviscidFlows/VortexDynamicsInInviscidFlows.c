#include <stdio.h>
#include <math.h>

#define PI 3.1415926535
#define N 41
#define GAM 1.0
#define r 0.6
#define a 0.5

float cacu_rm(float, float, float, float);      //Get Rm
float cacu_rn(float, float, float, float);		//Get Rn

float cacu_Um(float, float, float);             //Get Um
float cacu_Vm(float, float, float);             //Get Vm

float cacu_New_Um(float, int);					//Get New Um
float cacu_New_Vm(float, int);					//Get New Um

float cacu_um(float, float, float);             //Get um
float cacu_vm(float, float, float);             //Get vm

float cacu_sigmaud(float, float, float);             //Get ud
float cacu_sigmavd(float, float, float);             //Get vd

float cacu_qm(float, float);                    //Get qm

float cacu_rhs_m(float, float, float, float);   //Get rhsm

float cacu_rhs_p(float, float);                 //Get rhsp

float init_x(float, float);                     //initialization  X of the circle
float init_y(float, float);                     //initialization  Y of the circle
float init_cp(float, float);                    //initialization  Control points of the circle

float ang_to_rad(float);                        //trans degree to radian

int cacu_nsubs(float,float);                  //Get number of n subs

float init_tempXY(float, float, float);       //Get xn and yn for each edge

float cacu_beta(float,float,float,float,float,float);

int main(void)
{
	float x[N], y[N], tempXi[N], tempYi[N], tempXj[N], tempYj[N];
	float xm[N], xj[N], ym[N], yj[N], x0, y0, xi[N], yi[N], xn[N], yn[N], tempsubs[N];
	float tempX[N], tempY[N];
	float Umj[N], Vmj[N], Um0[N], Vm0[N], Umi[N], Vmi[N];
	float NewUmj[N], NewVmj[N];

	float um[N], vm[N];
	float Newum[N], Newvm[N];
	float ud[N], vd[N];

	float rmj[N], rm0[N], rmi[N], rnj[N]={0,};
	float qm[N];
	float Newqm[N];

	float ir;
	float deltaSn[N], deltaSp[N], deltaSm[N];
	float rhsm=0, rhsM[N], rhsP=0;
	float ceta, deltaBeta[N], cosbeta[N], sinbeta[N];
	float radian, angle;
	float sigmaU[N]={0,}, sigmaV[N]={0,};
	float sigmaud=0, sigmavd=0;

	int m, i, j, nsubs[N];
	
	FILE *test=fopen("test.dat","wt");
	FILE *test1=fopen("test1.dat","wt");
	FILE *test2=fopen("test2.dat","wt");

	ir = (float)pow(a,2) / r;
	x0=0.0, y0=0.0;

	for(i = 0, angle = 0; i <= N-1, angle <= 360; i++, angle += (360/(N-1)))
	{
		x[i] = init_x(ang_to_rad(angle),a);
		y[i] = init_y(ang_to_rad(angle),a);

		tempXi[i] = init_x(ang_to_rad(angle),ir);
		tempYi[i] = init_y(ang_to_rad(angle),ir);

		tempXj[i] = init_x(ang_to_rad(angle),r);
		tempYj[i] = init_y(ang_to_rad(angle),r);

		printf( "%0.5f  %0.5f\n", x[i], y[i] );

		fprintf(test, "%0.5f  %0.5f\n", x[i], y[i] );
	}

	for(i = 0; i < N - 1; i++)
	{
		xm[i] = init_cp(x[i+1],x[i]);
		ym[i] = init_cp(y[i+1],y[i]);

//		printf( "%0.5f  %0.5f\n", x[i], y[i] );

		xi[i] = init_cp(tempXi[i+1],tempXi[i]);
		yi[i] = init_cp(tempYi[i+1],tempYi[i]);

		xj[i] = init_cp(tempXj[i+1],tempXj[i]);
		yj[i] = init_cp(tempYj[i+1],tempYj[i]);

		deltaBeta[i] = cacu_beta(x[i+1],y[i+1],x0,y0,x[i],y[i]);
		deltaSn[i] = r - a;

		cosbeta[i] = (float)cos(ang_to_rad(deltaBeta[i]));
		sinbeta[i] = (float)sin(ang_to_rad(deltaBeta[i]));

		deltaSp[i] = (float)sqrt((float)pow((x[i+1] - x[i]),2) + (float)pow((y[i+1] - y[i]),2));
		deltaSm[i] = deltaSp[i];
		
		fprintf(test1, "%0.5f  %0.5f\n", xm[i], ym[i] );
	
//		printf( "%d  %0.5f  %0.5f  %0.5f\n", i+1, cosbeta[i], sinbeta[i], ang_to_rad(deltaBeta[i]) );
	}
	
	
//	xj[7] = r * cos(ang_to_rad(90));
//	yj[7] = r * sin(ang_to_rad(90));

//	xi[7] = ir * cos(ang_to_rad(90));
//	yi[7] = ir * sin(ang_to_rad(90));

	for(i = 0; i < N - 1; i++)
	{
		rm0[i] = cacu_rm(xm[i], x0, ym[i], y0);
		rmi[i] = cacu_rm(xm[i], xi[7], ym[i], yi[7]);
		rmj[i] = cacu_rm(xm[i], xj[7], ym[i], yj[7]);

		nsubs[i] = cacu_nsubs(rmj[i], deltaSm[i]);
		
		tempsubs[i] = (1/(float)(nsubs[i]));

	//	printf("\ndi %d duan:\n", i+1);

		for(j=0; j<nsubs[i]; j++)
		{	
			if(nsubs[i]==1)
			{
				tempX[j] = x[i];
				tempY[j] = y[i];
				tempX[j+1] = x[i+1];
				tempY[j+1] = y[i+1];
			}
			else if( (nsubs[i]!=1) && (j!=(nsubs[i]-1)) )
			{
				tempX[0] = x[i];
				tempY[0] = y[i];		
				tempX[j+1] = init_tempXY(x[i+1], x[i], tempsubs[i]);
				tempY[j+1] = init_tempXY(y[i+1], y[i], tempsubs[i]);
				tempsubs[i] += (1/(float)(nsubs[i]));
				tempX[j+2] = x[i+1];
				tempY[j+2] = y[i+1];
			}
			else
				break;
		}
		for(j=0; j<nsubs[i]; j++)
		{
			xn[j] = init_cp(tempX[j+1], tempX[j]);
			yn[j] = init_cp(tempY[j+1], tempY[j]);
			rnj[j] = cacu_rn(xn[j], xj[7], yn[j], yj[7]);
			sigmaU[i] += ((yn[j]-yj[7])/(float)pow(rnj[j],2));
			sigmaV[i] += ((xn[j]-xj[7])/(float)pow(rmj[j],2));
			fprintf(test2, "%0.5f  %0.5f\n", xn[j], yn[j] );
		}
	//	printf( "%0.5f  %0.5f\n", sigmaU[i], sigmaV[i] );
	//	printf( "%0.5f  %0.5f\n", rmj[i], rnj[i] );

		NewUmj[i] = cacu_New_Um(sigmaU[i], nsubs[i]);
		NewVmj[i] = cacu_New_Vm(sigmaV[i], nsubs[i]);

		Um0[i] = cacu_Um(ym[i], y0, rm0[i]);
		Umi[i] = cacu_Um(ym[i], yi[7], rmi[i]);
		Umj[i] = cacu_Um(ym[i], yj[7], rmj[i]);

//		printf( "%0.5f  %0.5f\n", NewUmj[i], Umj[i] );

		Vm0[i] = cacu_Vm(xm[i], x0, rm0[i]);
		Vmi[i] = cacu_Vm(xm[i], xi[7], rmi[i]);
		Vmj[i] = cacu_Vm(xm[i], xj[7], rmj[i]);

		rhsM[i] = cacu_rhs_m(Umj[i], cosbeta[i], Vmj[i], sinbeta[i]);

		um[i] = cacu_um(Um0[i], Umi[i], Umj[i]);
		vm[i] = cacu_vm(Vm0[i], Vmi[i], Vmj[i]);
		
		sigmaud += cacu_sigmaud(GAM, deltaSm[i], Umj[i]);
		sigmavd += cacu_sigmavd(GAM, deltaSm[i], Vmj[i]);

//		sigmaud += cacu_sigmaud(GAM, deltaSm[i], NewUmj[i]);
//		sigmavd += cacu_sigmavd(GAM, deltaSm[i], NewVmj[i]);

		qm[i] = cacu_qm(um[i], vm[i]);

		Newum[i] = cacu_um(Um0[i], Umi[i], NewUmj[i]);
		Newvm[i] = cacu_vm(Vm0[i], Vmi[i], NewVmj[i]);

		Newqm[i] = cacu_qm(Newum[i], Newvm[i]);

//	printf("%d um = %0.4f vm = %0.4f qm = %0.4f rhsm = %0.4f rmj = %0.4f nsubs = %d\n",
//		i+1, um[i], vm[i], qm[i], rhsM[i], rmj[i], nsubs[i]);
//	printf("%d um = %0.4f vm = %0.4f qm = %0.4f rhsm = %0.4f rmj = %0.4f nsubs = %d\n",
//		i+1, Newum[i], Newvm[i], Newqm[i], rhsM[i], rmj[i], nsubs[i]);
//	printf("%d   um = %0.4f   vm = %0.4f   qm = %0.4f\n", i+1, um[i], vm[i], qm[i]);
//		printf("%d  %0.4f  %0.4f  %0.4f  %0.4f\n", i+1, Um0[i],Umi[i],Umj[i],Um0[i]-Umi[i]+Umj[i]);
	}

	ud[0] = (-1) * sigmaud;
	vd[0] = (-1) * sigmavd;
	
//	printf("%0.4f  %0.4f\n", ud[0], vd[0] );

	for(i = 0; i < N - 1; i++)
	{
		rhsm += (rhsM[i] * deltaSn[i]);
	}

	rhsP = cacu_rhs_p(deltaSp[i], rhsm);

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

float cacu_New_Um(float PLSTOTAL, int NSUB_)
{
	float temp_new_Um_, temp_tot_;
	int temp_nsubs_;
	temp_tot_=PLSTOTAL; temp_nsubs_=NSUB_;
	temp_new_Um_ = (1/(2*PI*temp_nsubs_))*PLSTOTAL;
	return temp_new_Um_;
	temp_tot_=0; temp_nsubs_=0;
}

float cacu_New_Vm(float PLSTOTAL, int NSUB_)
{
	float temp_new_Vm_, temp_tot_;
	int temp_nsubs_;
	temp_tot_=PLSTOTAL; temp_nsubs_=NSUB_;
	temp_new_Vm_ = (-1/(2*PI*temp_nsubs_))*PLSTOTAL;
	return temp_new_Vm_;
	temp_tot_=0; temp_nsubs_=0;
}

float cacu_um(float M0, float Mi, float Mj)
{
	float temp_um_, temp_m0_, temp_mi_, temp_mj_;
	temp_m0_=M0; temp_mi_=Mi; temp_mj_=Mj;

	temp_um_=(temp_m0_ - temp_mi_ + temp_mj_) * GAM;
	return temp_um_;
	temp_um_=0; temp_m0_=0; temp_mi_=0; temp_mj_=0;
}

float cacu_vm(float M0, float Mi, float Mj)
{
	float temp_vm_, temp_m0_, temp_mi_, temp_mj_;
	temp_m0_=M0; temp_mi_=Mi; temp_mj_=Mj;

	temp_vm_=(temp_m0_ - temp_mi_ + temp_mj_) * GAM;
	return temp_vm_;
	temp_vm_=0; temp_m0_=0; temp_mi_=0; temp_mj_=0;
}

float cacu_sigmaud(float SM, float deltaSM, float Umj)
{
	float temp_sigmaud_, temp_sm_, temp_deltasm_, temp_umj_;
	temp_sm_=SM; temp_deltasm_=deltaSM; temp_umj_=Umj;
	temp_sigmaud_ = temp_sm_*temp_deltasm_*temp_umj_;
	return temp_sigmaud_;
	temp_sigmaud_=0; temp_sm_=0; temp_deltasm_=0; temp_umj_=0;
}	

float cacu_sigmavd(float SM, float deltaSM, float Vmj)
{
	float temp_sigmavd_, temp_sm_, temp_deltasm_, temp_vmj_;
	temp_sm_=SM; temp_deltasm_=deltaSM; temp_vmj_=Vmj;
	temp_sigmavd_ = temp_sm_*temp_deltasm_*temp_vmj_;
	return temp_sigmavd_;
	temp_sigmavd_=0; temp_sm_=0; temp_deltasm_=0; temp_vmj_=0;
}

float cacu_rm(float XM, float X_, float YM, float Y_)
{
	float temp_rm_, temp_xm_, temp_x_, temp_ym_, temp_y_;
	temp_xm_=XM; temp_x_=X_; temp_ym_=YM; temp_y_=Y_;
	temp_rm_ = (float)sqrt( (float)pow((temp_xm_-temp_x_),2) + (float)pow((temp_ym_-temp_y_),2) );
	return temp_rm_;
	temp_xm_=0; temp_x_=0;temp_ym_=0;temp_y_=0;temp_rm_=0;
}

float cacu_rn(float XN, float XJ, float YN, float YJ)
{
	float temp_rn_, temp_xn_, temp_xj_, temp_yn_, temp_yj_;
	temp_xn_=XN; temp_xj_=XJ; temp_yn_=YN; temp_yj_=YJ;
	temp_rn_ = (float)sqrt( (float)pow((temp_xn_-temp_xj_),2) + (float)pow((temp_yn_-temp_yj_),2) );
	return temp_rn_;
	temp_xn_=0; temp_xj_=0;temp_yn_=0;temp_yj_=0;temp_rn_=0;
}

float cacu_qm(float UM, float VM)
{
	float temp_q_, temp_um_, temp_vm_;
	temp_um_=UM; temp_vm_=VM;
	temp_q_ = (float)sqrt( (float)pow(temp_um_,2) + (float)pow(temp_vm_,2) );
	return temp_q_;
	temp_q_=0; temp_um_=0; temp_vm_=0;
}

float cacu_rhs_m(float Umj, float cosBm, float Vmj, float sinBm)
{
	float temp_rhsm_, temp_umj_, temp_cosbm_, temp_vmj_, temp_sinbm_, temp_gam_ = GAM;
	temp_umj_=Umj; temp_cosbm_=cosBm; temp_vmj_=Vmj; temp_sinbm_=sinBm;
	temp_rhsm_ = (-1) * temp_gam_ * (temp_umj_*temp_cosbm_ +  temp_vmj_*temp_sinbm_);
	return temp_rhsm_;
	temp_rhsm_=0; temp_umj_=0; temp_cosbm_=0; temp_vmj_=0; temp_sinbm_=0;
}

float cacu_rhs_p(float deltaSP, float rhsnPLUS)
{
	float temp_rhsp_, temp_deltasp_, temp_rhsn_;
	temp_deltasp_=deltaSP; temp_rhsn_=rhsnPLUS;
	temp_rhsp_ = (-1) * (1 / (temp_deltasp_)) * temp_rhsn_;
	return temp_rhsp_;
	temp_deltasp_=0; temp_rhsn_=0; temp_rhsp_=0;
}

float init_x(float cta, float R)
{
	float temp_x;
	temp_x=R*(float)cos(cta);
	return temp_x;
}

float init_y(float cta, float R)
{
	float temp_y;
	temp_y=R*(float)sin(cta);
	return temp_y;
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

int cacu_nsubs(float RMJ, float DELTASM)
{
	float temp_nsubs_, temp_rmj_, temp_deltasm_;
	int n;
	temp_rmj_=RMJ; temp_deltasm_=DELTASM;
	temp_nsubs_=1+((2*temp_deltasm_)/temp_rmj_);
	n=(int)temp_nsubs_;
	return n;
	temp_nsubs_=0;temp_rmj_=0;temp_deltasm_=0;n=0;	
}

float init_tempXY(float V_1, float V, float SUBS)
{
	float temp_xy, temp_subs_, temp_v1_, temp_v_;
	temp_subs_ = SUBS; temp_v1_ = V_1; temp_v_ = V;
	if(temp_subs_ != 1)
	{
		temp_xy = temp_subs_ * (temp_v1_ - temp_v_) + temp_v_;
	}
	return temp_xy;
	temp_xy=0; temp_subs_=0; temp_v1_ = 0; temp_v_ = 0;
}

float cacu_beta(float X1, float Y1, float X0, float Y0, float X2, float Y2)
{
	float temp_x1_, temp_y1_, temp_x0_, temp_y0_, temp_x2_, temp_y2_;
	float temp_rad_, temp_beta_;
	temp_x1_=X1; temp_y1_=Y1; temp_x0_=X0; temp_y0_=Y0; temp_x2_=X2; temp_y2_=Y2;
	temp_rad_ = (float)fabs((float)atan((temp_y0_-temp_y1_)/(temp_x0_-temp_x1_))-(float)atan((temp_y2_-temp_y0_)/(temp_x2_-temp_x0_)));
	temp_beta_ = temp_rad_ * (180 / PI);
	if(((temp_x1_<0) && (temp_x2_>0))||((temp_x1_>0) && (temp_x2_<0)))
		return 180 - temp_beta_;
	else
		return temp_beta_;
	temp_x1_=0; temp_y1_=0; temp_x0_=0; temp_y0_=0; temp_x2_=0; temp_y2_=0;
	temp_rad_=0, temp_beta_=0;
}
