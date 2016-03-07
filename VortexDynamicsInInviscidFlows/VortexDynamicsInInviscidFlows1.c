#include <stdio.h>
#include <math.h>

#define PI 3.1415926535
#define GAM 1.0
#define N 2

float cacu_Umn(float, float, float, float);
float cacu_Vmn(float, float, float, float);
float cacu_u(float);
float cacu_v(float);
float xb(float, float, float);
float yb(float, float, float);
float xd(float, float, float, float);
float yd(float, float, float, float);

int main(void)
{
	float U[N], V[N];
	float u[N], v[N];
	float x[N], y[N], xb[N], yb[N], xd[N], yd[N];

	float sigmaU[N], sigmaV[N];

	float deltaT;

	FILE *test1 = fopen("test1.dat","wt");
	FILE *test2 = fopen("test2.dat","wt");

	x[0] = 0.0; y[0] = 0.5;
	x[1] = 0.0; y[1] = -0.5;


	return 0;
}

float cacu_Umn(float Xm, float Ym, float Xn, float Yn)
{
	float temp_umn_, temp_xm_, temp_ym_, temp_xn_, temp_yn_;
	temp_xm_=Xm; temp_ym_=Ym; temp_xn_=Xn; temp_yn_=Yn;

	temp_umn_ = (1/(2*PI)) * ((temp_ym_-temp_yn_) / (pow((temp_xm_-temp_xn_),2) + pow((temp_ym_-temp_yn_),2)));

	return temp_umn_;
	temp_xm_=0; temp_ym_=0; temp_xn_=0; temp_yn_=0;
}

float cacu_Vmn(float Xm, float Ym, float Xn, float Yn)
{
	float temp_vmn_, temp_xm_, temp_ym_, temp_xn_, temp_yn_;
	temp_xm_=Xm; temp_ym_=Ym; temp_xn_=Xn; temp_yn_=Yn;

	temp_vmn_ = (1/(2*PI)) * ((temp_xm_-temp_xn_) / (pow((temp_xm_-temp_xn_),2) + pow((temp_ym_-temp_yn_),2)));

	return temp_vmn_;
	temp_vmn_=0; temp_xm_=0; temp_ym_=0; temp_xn_=0; temp_yn_=0;
}

float cacu_u(float Umn)
{
	float temp_u_, temp_umn_;
	temp_umn_=Umn;

	temp_u_ = GAM * temp_umn_;

	return temp_u_;
	temp_u_=0; temp_umn_=0;
}

float cacu_v(float Vmn)
{
	float temp_v_, temp_vmn_;
	temp_vmn_=Vmn;

	temp_v_ = GAM * temp_vmn_;

	return temp_v_;
	temp_v_=0; temp_vmn_=0;
}

float xb(float Xa, float u, float dT)
{
	float temp_xb_, temp_xa_, temp_u_, temp_dt_;
	temp_xa_=Xa; temp_u_=u; temp_dt_=dT;

	temp_xb_ = temp_xa_ + temp_u_ * temp_dt_;

	return temp_xb_;
	temp_xb_=0; temp_xa_=0; temp_u_=0; temp_dt_=0;
}

float yb(float Ya, float v, float dT)
{
	float temp_yb_, temp_ya_, temp_v_, temp_dt_;
	temp_ya_=Ya; temp_v_=v; temp_dt_=dT;

	temp_yb_ = temp_ya_ + temp_v_ * temp_dt_;

	return temp_yb_;
	temp_yb_=0; temp_ya_=0; temp_v_=0; temp_dt_=0;
}

float xd(float Xa, float ua, float ub, float dT)
{
	float temp_xd_, temp_xa_, temp_ua_, temp_ub_, temp_dt_;
	temp_xa_=Xa; temp_ua_=ua; temp_ub_=ub; temp_dt_=dT;

	temp_xd_ = temp_xa_ + 0.5 * (temp_ua_ + temp_ub_) * temp_dt_;

	return temp_xd_;
	temp_xd_=0; temp_xa_=0; temp_ua_=0; temp_ub_=0; temp_dt_=0;
}

float yd(float Ya, float va, float vb, float dT)
{
	float temp_yd_, temp_ya_, temp_va_, temp_vb_, temp_dt_;
	temp_ya_=Ya; temp_va_=va; temp_vb_=vb; temp_dt_=dT;

	temp_yd_ = temp_ya_ + 0.5 * (temp_va_ + temp_vb_) * temp_dt_;

	return temp_yd_;
	temp_yd_=0; temp_ya_=0; temp_va_=0; temp_vb_=0; temp_dt_=0;
}
