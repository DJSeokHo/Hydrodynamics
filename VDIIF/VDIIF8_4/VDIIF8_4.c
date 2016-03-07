#include <stdio.h>
#include <math.h>

#define PI 3.1415926535
#define N 41
#define GAM 1.0
#define a 0.5
#define ratX 1.0
#define ratY 1.0

float cacu_beta(float,float,float,float,float,float);	//Get angle beta

float cacu_rhs_m(float, float, float, float);   //Get rhsm
float cacu_rhs_p(float, float);                 //Get rhsp

float cacu_rm(float, float, float, float);      //Get Rm
float cacu_rn(float, float, float, float);		//Get Rn

float cacu_Um(float, float, float);             //Get Um
float cacu_Vm(float, float, float);             //Get Vm

float cacu_New_Um(float, int);					//Get New Um
float cacu_New_Vm(float, int);					//Get New Um

float cacu_um(float, float, float);             //Get um
float cacu_vm(float, float, float);             //Get vm

float cacu_qm(float, float);                    //Get qm

float cacu_t(float, float);						//Get time

float cacu_tempr(float, float, float, float);	//Get temp radiu

int cacu_nsubs(float,float);					//Get number of sub-elements

float init_x(float, float);                     //initialization  X of the circle
float init_y(float, float);                     //initialization  Y of the circle
float init_cp(float, float);                    //initialization  Control points of the circle

float ang_to_rad(float);                        //trans degree to radian
float rad_to_ang(float);
float init_tempXY(float, float, float);       //Get xn and yn for each edge

float cacu_xb(float, float, float);
float cacu_yb(float, float, float);

float cacu_xd(float, float, float, float);
float cacu_yd(float, float, float, float);

float init_r(float, float, float, float);

main()
{
	float x[N], y[N], tempXi[N], tempYi[N], tempXj[N], tempYj[N];
	float xm[N], ym[N], xja, yja, xjb, yjb, xjd, yjd, x0, y0, xi, yi, xn[N], yn[N], tempsubs[N];

	float tempX[N], tempY[N];
	float Umj[N], Vmj[N], Um0[N], Vm0[N], Umi[N], Vmi[N];
	float NewUmj[N], NewVmj[N];

	float um[N]={0,}, vm[N]={0,};
	float Newum[N], Newvm[N];

	float ub = 0, vb = 0;
	float ud = 0, vd = 0;

	float rmj[N], rm0[N], rmi[N], rnj[N]={0,};
	float qm[N];
	float Newqm[N];

	float realR, realA, ir;
	float deltaSn[N], deltaSm[N];
	float NewrhsM[N], rhsM[N], rhsN[N], rhsm=0, rhsP=0;
	float ceta, deltaBeta[N], cosbeta[N], sinbeta[N];
	float radian, angle;
	float sigmaU[N]={0,}, sigmaV[N]={0,};
	float sigmaud=0, sigmavd=0;
	
	float X[N]={0,}, Y[N]={0,};
	
	float ElementLength[N]={0,}, deltaSp=0.1045;
	
	float RatioOfGap;
	
	float umtest=0;

	float t, deltaT;

	float gamSM[N];

	float tempR, tempx, tempy;

	float setIcosceta=0, setIsinceta=0, tempAngleSin=0, tempAngleCos=0, tempAngle=0;

	int m, i, j, nsubs[N], time;

	FILE *circle=fopen("circle.dat","wt");
	FILE *pivotalP=fopen("pivotalP.dat","wt");
	FILE *inversionP=fopen("inversionP.dat","wt");
	FILE *sub=fopen("sub.dat","wt");
	FILE *test=fopen("test.dat","wt");
	FILE *driftPath=fopen("driftPath.dat","wt");

	for(i = 0, angle = 0; i <= N-1, angle <= 360; i++, angle += (360/(N-1)))
	{
		x[i] = init_x(ang_to_rad(angle),a * ratX);
		y[i] = init_y(ang_to_rad(angle),a * ratY);
		fprintf(circle, "%0.4f  %0.4f\n", x[i], y[i] );
	}

	for(i = 0; i < N - 1; i++)
	{
		xm[i] = init_cp(x[i+1],x[i]);
		ym[i] = init_cp(y[i+1],y[i]);

		fprintf(pivotalP, "%0.4f  %0.4f\n", xm[i], ym[i] );
	}

	for(i = 0; i < N - 1; i++)	//element's deltaS
	{
		deltaBeta[i] = cacu_beta(x[i+1],y[i+1],x0,y0,x[i],y[i]);	//delta angle beta
		cosbeta[i] = (float)cos(ang_to_rad(deltaBeta[i]));
		sinbeta[i] = (float)sin(ang_to_rad(deltaBeta[i]));

		ElementLength[i] = (float)sqrt((float)pow((x[i+1] - x[i]),2) + (float)pow((y[i+1] - y[i]),2));
		deltaSm[i] = ElementLength[i];
	}

//	xm[7]=0.0; ym[7]=0.5;				//use this can get the right value as report shown

	xja=0.0; yja=0.6;
	x0=0.0; y0=0.0;
	xi=0.0; yi=0.0;

	deltaT = 0.12506;
	fprintf(driftPath,"%0.4f %0.4f\n", xja,yja);

	for(time=0; time<40; time++)
	{
		realR = init_r(x0, y0, xja, yja);			//(X0, Y0, X1, Y1)
				
		setIcosceta = xja/realR;					//This is ceta of (xj,yj)
		setIsinceta = yja/realR;					//This is ceta of (xj,yj)
			
		tempAngleCos = (float)acos(setIcosceta);		//Get angle of (xj,yj)
		tempAngleSin = (float)asin(setIsinceta);		//Get angle of (xj,yj)

		tempx = init_x(tempAngleCos, a * ratX);		//Calculate (x,y) of circle
		tempy = init_y(tempAngleSin, a * ratY);		//Calculate (x,y) of circle

		realA = init_r(x0, y0, tempx, tempy);		//Calculate radiu of circle

		ir = (float)pow(realA,2) / realR;			//Calculate radiu of (xi,yi)

		xi = init_x(tempAngleCos, ir);					//Calculate (x,y) of point i
		yi = init_y(tempAngleSin, ir);					//Calculate (x,y) of point i

	//	fprintf(inversionP,"%0.4f %0.4f\n", xi, yi);
		
		for(i = 0; i < N - 1; i++)
		{
			rm0[i] = cacu_rm(xm[i], x0, ym[i], y0);
			rmi[i] = cacu_rm(xm[i], xi, ym[i], yi);
			rmj[i] = cacu_rm(xm[i], xja, ym[i], yja);

			Um0[i] = cacu_Um(ym[i], y0, rm0[i]);
			Umi[i] = cacu_Um(ym[i], yi, rmi[i]);
			Umj[i] = cacu_Um(ym[i], yja, rmj[i]);

			Vm0[i] = cacu_Vm(xm[i], x0, rm0[i]);
			Vmi[i] = cacu_Vm(xm[i], xi, rmi[i]);
			Vmj[i] = cacu_Vm(xm[i], xja, rmj[i]);
			
	//		rhsM[i] = cacu_rhs_m(Umj[i], cosbeta[i], Vmj[i], sinbeta[i]);
	//		rhsN[i] = rhsM[i];
	//		rhsm += rhsN[i] * deltaSn[i];

			um[i] = cacu_um(Um0[i], Umi[i], Umj[i]);
			vm[i] = cacu_vm(Vm0[i], Vmi[i], Vmj[i]);
			qm[i] = cacu_qm(um[i], vm[i]);

			gamSM[i] = qm[i];

			ub += (gamSM[i]*deltaSm[i]*Umj[i]);
			vb += (gamSM[i]*deltaSm[i]*Vmj[i]);

		//	printf("%0.4f  %0.4f\n", ub, vb);

	/*		nsubs[i] = cacu_nsubs(rmj[i], ElementLength[i]);
			tempsubs[i] = (1/(float)(nsubs[i]));

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
				rnj[j] = cacu_rn(xn[j], xj, yn[j], yj);
				sigmaU[i] += ((yn[j]-yj)/(float)pow(rnj[j],2));
				sigmaV[i] += ((xn[j]-xj)/(float)pow(rnj[j],2));
				fprintf(sub, "%0.5f  %0.5f\n", xn[j], yn[j] );
			}

			NewUmj[i] = cacu_New_Um(sigmaU[i], nsubs[i]);
			NewVmj[i] = cacu_New_Vm(sigmaV[i], nsubs[i]);

			NewrhsM[i] = cacu_rhs_m(NewUmj[i], cosbeta[i], NewVmj[i], sinbeta[i]);

			printf("%d -- %0.4f\n", i+1, rmj[i]);
			printf("%d -- %0.4f, %0.4f, %0.4f, %0.4f\n", i+1, rmj[i], um[i], vm[i], qm[i]);
	*/
		}

	//	ud = (-1)*ud;
	//	vd = (-1)*vd;
	//	printf("%0.4f  %0.4f\n", ud, vd);
		
		xjb = xja + ub * deltaT;
		yjb = yja + vb * deltaT;

//		fprintf(inversionP,"%0.4f %0.4f\n", xi, yi);
//		fprintf(test,"%0.4f %0.4f\n", xja,yja);
//		fprintf(test,"%0.4f %0.4f\n", xjb,yjb);
		realR = 0; setIcosceta=0; setIsinceta=0; 
		tempx = 0; tempy = 0; realA = 0; ir = 0;
		xi = 0; yi = 0;

		realR = init_r(x0, y0, xjb, yjb);			//(X0, Y0, X1, Y1)
				
		setIcosceta = xjb/realR;					//This is ceta of (xj,yj)
		setIsinceta = yjb/realR;					//This is ceta of (xj,yj)
			
		tempAngleCos = (float)acos(setIcosceta);		//Get angle of (xj,yj)
		tempAngleSin = (float)asin(setIsinceta);		//Get angle of (xj,yj)

		tempx = init_x(tempAngleCos, a * ratX);		//Calculate (x,y) of circle
		tempy = init_y(tempAngleSin, a * ratY);		//Calculate (x,y) of circle

		realA = init_r(x0, y0, tempx, tempy);		//Calculate radiu of circle

		ir = (float)pow(realA,2) / realR;			//Calculate radiu of (xi,yi)

		xi = init_x(tempAngleCos, ir);					//Calculate (x,y) of point i
		yi = init_y(tempAngleSin, ir);	
//		ub=0; vb=0;
		
//		xja=xjb;yja=yjb;

		for(i = 0; i < N - 1; i++)
		{
			rm0[i] = cacu_rm(xm[i], x0, ym[i], y0);
			rmi[i] = cacu_rm(xm[i], xi, ym[i], yi);
			rmj[i] = cacu_rm(xm[i], xjb, ym[i], yjb);

			Um0[i] = cacu_Um(ym[i], y0, rm0[i]);
			Umi[i] = cacu_Um(ym[i], yi, rmi[i]);
			Umj[i] = cacu_Um(ym[i], yjb, rmj[i]);

			Vm0[i] = cacu_Vm(xm[i], x0, rm0[i]);
			Vmi[i] = cacu_Vm(xm[i], xi, rmi[i]);
			Vmj[i] = cacu_Vm(xm[i], xjb, rmj[i]);
			
	//		rhsM[i] = cacu_rhs_m(Umj[i], cosbeta[i], Vmj[i], sinbeta[i]);
	//		rhsN[i] = rhsM[i];
	//		rhsm += rhsN[i] * deltaSn[i];

			um[i] = cacu_um(Um0[i], Umi[i], Umj[i]);
			vm[i] = cacu_vm(Vm0[i], Vmi[i], Vmj[i]);
			qm[i] = cacu_qm(um[i], vm[i]);

			gamSM[i] = qm[i];

			ud += (gamSM[i]*deltaSm[i]*Umj[i]);
			vd += (gamSM[i]*deltaSm[i]*Vmj[i]);

	//		printf("%0.4f  %0.4f\n", ud, vd);

	/*		nsubs[i] = cacu_nsubs(rmj[i], ElementLength[i]);
			tempsubs[i] = (1/(float)(nsubs[i]));

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
				rnj[j] = cacu_rn(xn[j], xj, yn[j], yj);
				sigmaU[i] += ((yn[j]-yj)/(float)pow(rnj[j],2));
				sigmaV[i] += ((xn[j]-xj)/(float)pow(rnj[j],2));
				fprintf(sub, "%0.5f  %0.5f\n", xn[j], yn[j] );
			}

			NewUmj[i] = cacu_New_Um(sigmaU[i], nsubs[i]);
			NewVmj[i] = cacu_New_Vm(sigmaV[i], nsubs[i]);

			NewrhsM[i] = cacu_rhs_m(NewUmj[i], cosbeta[i], NewVmj[i], sinbeta[i]);

			printf("%d -- %0.4f\n", i+1, rmj[i]);
			printf("%d -- %0.4f, %0.4f, %0.4f, %0.4f\n", i+1, rmj[i], um[i], vm[i], qm[i]);
	*/

		}

		xjd = xja + 0.5*(ub+ud) * deltaT;
		yjd = yja + 0.5*(vb+vd) * deltaT;
		
	//	fprintf(test,"%0.4f %0.4f\n", xjb,yjb);
		fprintf(driftPath,"%0.4f %0.4f\n", xjd,yjd);
		
		xja=xjd;yja=yjd;
		ub=0;vb=0;ud=0;vd=0;
		realR = 0; setIcosceta=0; setIsinceta=0; 
		tempx = 0; tempy = 0; realA = 0; ir = 0;
		xi = 0; yi = 0;
	}
//	printf("%0.4f  %0.4f\n", xja, yja);
//	printf("%0.4f  %0.4f\n", xjb, yjb);

//	rhsP=((-1)/deltaSp*rhsm);

	return 0;
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

float cacu_t(float R, float A)
{
	float temp_t_, temp_qd_, temp_r_, temp_a_;
	temp_r_=R; temp_a_=A;
	temp_qd_=(GAM/(2*PI*temp_r_))*((float)pow(temp_a_,2)/((float)pow(temp_r_,2)-(float)pow(temp_a_,2)));
	temp_t_=(2*PI*temp_r_)/temp_qd_;
	return temp_t_;
	temp_t_=0; temp_qd_=0; temp_r_=0; temp_a_=0;
}

float cacu_tempr(float X0, float Y0, float XJ, float YJ)
{
	float temp_r_, temp_x0_, temp_y0_, temp_xj_, temp_yj_;
	temp_x0_=X0; temp_y0_=Y0; temp_xj_=XJ; temp_yj_=YJ;
	temp_r_=(float)sqrt((float)pow((temp_x0_-temp_xj_),2)+(float)pow((temp_y0_-temp_yj_),2));
	return temp_r_;
	temp_r_=0; temp_x0_=0; temp_y0_=0; temp_xj_=0; temp_yj_=0;
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

float init_r(float X0, float Y0, float X1, float Y1)
{
	float temp_r_, temp_x0_, temp_y0_, temp_x1_, temp_y1_;
	temp_x0_=X0; temp_y0_=Y0; temp_x1_=X1; temp_y1_=Y1;
	temp_r_=(float)sqrt((float)pow((temp_x0_-temp_x1_),2)+(float)pow((temp_y0_-temp_y1_),2));
	return temp_r_;
	temp_r_=0; temp_x0_=0; temp_y0_=0; temp_x1_=0; temp_y1_=0;
}

float rad_to_ang(float RAD)
{
	float temp_ang_, temp_rad_;
	temp_rad_=RAD;
	temp_ang_ = temp_rad_/PI*180;
	return temp_ang_;
	temp_rad_=0; temp_ang_=0;
}