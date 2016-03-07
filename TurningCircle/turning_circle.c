//***********************************************************************************//
//         The main dimensions of a general purpose cargo vessel are:                //
//             Length                      L = 140.00 m								 //
//             Breadth                     B = 18.20 m								 //
//             Draught                     d = 7.82 m								 //
//             Block coefficient           Cb = 0.673								 //
//			   Projected rudder area       AR = 14.60 m^2							 //
//             Initial speed               17.0 konts(1 kont = 0.5144 m/s)			 //
//			   rudder angle                20										 //
//             rudder speed                2.0/s									 //
//             t1                          20/2                                      //
//             t2                          460 s									 //
//             dPsi                        0.90/s                                    //
//             dropSpeed                   11.9 knots	                             //
//             driftAngle                  8.2										 //
//																					 //
//							Code--By--SeokHo(SHI HAO)								 //									
//																					 //
//***********************************************************************************//

#include <iostream>
#include <cmath>

#define PI 3.1415926535
#define L 140.00
#define G 70.00
#define B 18.20
#define d 7.82
#define Cb 0.673
#define AR 14.60
#define InitV 17.0
#define rudA 20
#define rudS 2.0
#define t1 20/2
#define t2 460
#define Psi 360
#define dPsi 0.90
#define dropspeed 11.9
#define driftAngle 8.2

using namespace std;

float ang_to_rad(float);
float rad_to_ang(float);

float Turning_circle_diameter(float, float);
float drift_velocity(float, float);
float Neutral_or_pivoting_point(float, float, float);
float dirft_velocity_near_rudder(float, float, float, float);
float dirft_angle_near_rudder_U(float, float);
float dirft_angle_near_rudder_Rc(float Lp, float Xr, float dropV, float dP);
float Effective_rudder_angle(float, float);
float Nomote_indices_K(float, float);
float Nomote_indices_T(float, float, float, float);
float Turning_c_d_for_delta(float, float, float);

int main(void)
{
	cout << "turning circle diameter is " << Turning_circle_diameter(dropspeed, dPsi) << endl;
	cout << "drift velocity is " << drift_velocity(dropspeed, driftAngle) << endl;
	cout << "Neutral or pivoting point is " << Neutral_or_pivoting_point(dropspeed, dPsi, driftAngle) << endl;
	cout << "dirft velocity near rudder is " << dirft_velocity_near_rudder(G, dropspeed, dPsi, driftAngle) << endl;
	cout << "dirft angle near rudder by U is " << dirft_angle_near_rudder_U(dirft_velocity_near_rudder(G, dropspeed, dPsi, driftAngle), dropspeed) << endl;
	cout << "dirft angle near rudder by Rc is " << dirft_angle_near_rudder_Rc(Neutral_or_pivoting_point(dropspeed, dPsi, driftAngle), G, dropspeed, dPsi) << endl;
	cout << "effective rudder angle is " << Effective_rudder_angle(rudA, dirft_angle_near_rudder_U(dirft_velocity_near_rudder(G, dropspeed, dPsi, driftAngle), dropspeed)) << endl;
	cout << "Nomoto coefficients K is " << Nomote_indices_K(rudA, dPsi) << endl;
	cout << "Nomoto coefficients T is " << Nomote_indices_T(Psi, dPsi, t1, t2) << endl;
	cout << "turning circle diameter by 30 is " << Turning_c_d_for_delta(Nomote_indices_K(rudA, dPsi), 30, dropspeed) << endl;
	return 0;
}

float Turning_circle_diameter(float dropV, float dP)
{
	float temp_Dc, temp_Rc, temp_dPsi, temp_dropV;
	temp_dPsi=dP; temp_dropV=dropV;
	temp_Rc = (temp_dropV*0.5144) / (temp_dPsi * (PI/180));
	temp_Dc = 2 * temp_Rc;
	return temp_Dc;
}

float drift_velocity(float dropV, float beta)
{
	float temp_v, temp_dropV, temp_beta;
	temp_dropV=dropV; temp_beta=beta;
	temp_v = (-1) * temp_dropV * 0.5144 * (float)sin(ang_to_rad(temp_beta));
	return temp_v;
}

float Neutral_or_pivoting_point(float dropV, float dP, float beta)
{
	float temp_Rc, temp_dPsi, temp_dropV, temp_Lp, temp_beta;
	temp_dPsi=dP; temp_dropV=dropV; temp_beta=beta;
	temp_Rc = (temp_dropV*0.5144) / (temp_dPsi * (PI/180));
	temp_Lp = temp_Rc * (float)tan(ang_to_rad(temp_beta));
	return temp_Lp;
}

float dirft_velocity_near_rudder(float Xr, float dropV, float dP, float beta)
{
	float temp_VR, temp_Xr, temp_dropV, temp_dPsi, temp_beta;
	temp_Xr=Xr; temp_dropV=dropV; temp_dPsi=dP; temp_beta=beta;
	temp_VR = ((-1) * temp_dropV * 0.5144 * (float)sin(ang_to_rad(temp_beta)))
		+ (-1) * (temp_Xr * temp_dPsi * (PI/180));
	return temp_VR;
}

float dirft_angle_near_rudder_U(float Vr, float dropV)
{
	float temp_betaR, temp_Vr, temp_dropV;
	temp_Vr=Vr; temp_dropV=dropV;
	temp_betaR = (float)atan(((-1)*temp_Vr)/(temp_dropV*0.5144));
	return rad_to_ang(temp_betaR);
}

float dirft_angle_near_rudder_Rc(float Lp, float Xr, float dropV, float dP)
{
	float temp_Rc, temp_betaR, temp_Lp, temp_Xr, temp_dropV, temp_dPsi;
	temp_Lp=Lp; temp_Xr=Xr; temp_dPsi=dP; temp_dropV=dropV;
	temp_Rc = (temp_dropV*0.5144) / (temp_dPsi * (PI/180));
	temp_betaR = (float)atan((temp_Lp + temp_Xr)/temp_Rc);
	return rad_to_ang(temp_betaR);
}

float Effective_rudder_angle(float delta, float betaR)
{
	float temp_deltaS, temp_delta, temp_betaR;
	temp_delta=delta; temp_betaR=betaR;
	temp_deltaS = temp_delta - temp_betaR;
	return temp_deltaS;
}

float Nomote_indices_K(float delta, float dP)
{
	float temp_K, temp_delta, temp_dPsi;
	temp_delta=delta; temp_dPsi=dP;
	temp_K = temp_dPsi / temp_delta;
	return temp_K;
}

float Nomote_indices_T(float P, float dP, float tone, float ttwo)
{
	float temp_T, temp_Psi, temp_dPsi, temp_t1, temp_t2;
	temp_Psi=P; temp_dPsi=dP; temp_t1=tone; temp_t2=ttwo;
	temp_T = (-1) * temp_Psi / temp_dPsi + temp_t2 - temp_t1 / 2;
	return temp_T;
}

float Turning_c_d_for_delta(float K, float delta, float dropV)
{
	float temp_Dc, temp_Psi, temp_K, temp_delta, temp_dropV;
	temp_K=K; temp_delta=delta; temp_dropV=dropV;
	temp_Psi = temp_K * temp_delta;
	temp_Dc = (2 * temp_dropV * 0.5144) / (temp_Psi * (PI/180));
	return temp_Dc;
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
