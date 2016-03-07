#include <stdio.h>
#include <math.h>

#define N 30															//number of points

float cacu_U(float, float, float, float);
float cacu_V(float, float, float, float);
float cacu_u(float);
float cacu_v(float);
float cacu_xb(float, float, float);
float cacu_yb(float, float, float);
float cacu_xd(float, float, float, float);
float cacu_yd(float, float, float, float);
float cacu_error(float, float, float, float, float);

void draw(void);

main(void)
{
	draw();
	return 0;
}

void draw(void)
{
	float Umn[N]={0,}, Vmn[N]={0,}, deltaT;

	float ua[N]={0,}, va[N]={0,};
	float ub[N]={0,}, vb[N]={0,};

	float xa[N]={0,}, ya[N]={0,};
	float xb[N]={0,}, yb[N]={0,};
	float xd[N]={0,}, yd[N]={0,};

	float timestep = 0;
	float iter = 0; 

	int i, j, time;

	FILE *sheet=fopen("sheet.dat","wt");
	FILE *vortex=fopen("vortex.dat","wt");

	deltaT = 0.001;
	
	for(i=0; i<N; i++)														//initialize line
	{
		xa[i] = (float)i * 0.033333 + 0.005;
		ya[i] = 0.0;
		if(i==((N/2)-1))
		{
		//	xa[i] = xa[i] + 0.01;
		//	ya[i] = ya[i] + 0.01;
		}
		if(i==(N/2))
		{
		//	xa[i] = xa[i] - 0.01;
		//	ya[i] = ya[i] - 0.01;
		}
		if(i==0)
		{
			fprintf(sheet, "\nTitle=¡°sheet¡±\n");
			fprintf(sheet, "Variables=¡°X¡±,¡°Y\n");
			fprintf(sheet, "Zone I = %d \n\n", (N/2) );
		}
		if(i==(N/2))
		{
			fprintf(sheet, "\nTitle=¡°sheet¡±\n");
			fprintf(sheet, "Variables=¡°X¡±,¡°Y\n");
			fprintf(sheet, "Zone I = %d \n\n", (N/2) );
		}
		fprintf(sheet, "%0.8f %0.8f\n", xa[i], ya[i]);
	}
	
	for(time=0;time<300;time++)
	{
		for(i=0; i<N; i++)
		{
			for(j=0;j<N;j++)
			{
				if(i!=j)
				{
					Umn[j] = cacu_U(xa[i], ya[i], xa[j], ya[j]);
					Vmn[j] = cacu_V(xa[i], ya[i], xa[j], ya[j]);
					ua[i] += cacu_u(Umn[j]);
					va[i] += cacu_v(Vmn[j]);
				}
			}

			xb[i] = cacu_xb(xa[i], ua[i], deltaT);
			yb[i] = cacu_yb(ya[i], va[i], deltaT);

			for(j=0;j<N;j++)
			{
				Umn[j]=0;  Vmn[j]=0; 
			}
		}

		for(i=0; i<N; i++)
		{
			for(j=0;j<N;j++)
			{
				if(i!=j)
				{
					Umn[j] = cacu_U(xb[i], yb[i], xb[j], yb[j]);
					Vmn[j] = cacu_V(xb[i], yb[i], xb[j], yb[j]);
					ub[i] += cacu_u(Umn[j]);
					vb[i] += cacu_v(Vmn[j]);
				
				}
			}

			xd[i] = cacu_xd(xa[i], ua[i], ub[i], deltaT);
			yd[i] = cacu_yd(ya[i], va[i], vb[i], deltaT);

			for(j=0;j<N;j++)
			{
				Umn[j]=0;  Vmn[j]=0; 
			}
		}

		for(i=0;i<N;i++)														//initialize line
		{
			xa[i] = xd[i];
			ya[i] = yd[i];

			ua[i] = 0;															//importent operator
			va[i] = 0;															//importent operator
			ub[i] = 0;															//importent operator
			vb[i] = 0;															//importent operator
		}
		fprintf( sheet, "\nTitle=¡°sheet¡±\n");
		fprintf( sheet, "Variables=¡°X¡±,¡°Y\n");
		fprintf( sheet, "Zone I = %d \n\n", N );
		for(i=0;i<N;i++)
		{
			fprintf(sheet, "%0.6f %0.6f\n", xd[i], yd[i]);
		}
	}
	for(i=0; i<N; i++)
	{
	//	fprintf(vortex, "%0.8f %0.8f\n", xd[i], yd[i]);
	}


}

float cacu_U(float Xm, float Ym, float Xn, float Yn)
{
	float GAM = 1.0, T = 1.0;
	float PI=3.1415926535;
	float temp_U, tempXm, tempYm, tempXn, tempYn;
	tempXm=Xm; tempYm=Ym; tempXn=Xn; tempYn=Yn;
	temp_U = (GAM/(2*T))*( (float)sinh(((2*PI)/T)*(tempYm-tempYn))/((float)cosh(((2*PI)/T)*(tempYm-tempYn))-(float)cos(((2*PI)/T)*(tempXm-tempXn))));
	return temp_U;
	temp_U=0; tempXm=0; tempYm=0; tempXn=0; tempYn=0;
}
float cacu_V(float Xm, float Ym, float Xn, float Yn)
{
	float GAM = 1.0, T = 1.0;
	float PI=3.1415926535;
	float temp_V, tempXm, tempYm, tempXn, tempYn;
	tempXm=Xm; tempYm=Ym; tempXn=Xn; tempYn=Yn;
	temp_V = (-1)*(GAM/(2*T))*( (float)sin(((2*PI)/T)*(tempXm-tempXn))/((float)cosh(((2*PI)/T)*(tempYm-tempYn))-(float)cos(((2*PI)/T)*(tempXm-tempXn))));
	return temp_V;
	temp_V=0; tempXm=0; tempYm=0; tempXn=0; tempYn=0;
}
float cacu_u(float U)
{
	float GAM=1.0;
	float PI=3.1415926535;
	float temp_u, tempU;
	tempU=U;
	temp_u=GAM * tempU;
	return temp_u;
	temp_u=0; tempU=0;
}
float cacu_v(float V)
{
	float GAM=1.0;
	float PI=3.1415926535;
	float temp_v, tempV;
	tempV=V;
	temp_v=GAM * tempV;
	return temp_v;
	temp_v=0; tempV=0;
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