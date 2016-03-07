#include <stdio.h>
#include <math.h>
#define N 151
#define GAMm 1.0
#define GAMn 1.0

float cacu_U(float, float, float, float);
float cacu_V(float, float, float, float);
float cacu_u(float, float);
float cacu_v(float, float);
float cacu_xb(float, float, float);
float cacu_yb(float, float, float);
float cacu_xd(float, float, float, float);
float cacu_yd(float, float, float, float);
float cacu_error(float, float, float, float, float);

void drawing(void);

main()
{
	drawing();
	return 0;
}

void drawing(void)
{
	FILE *m=fopen("m.dat","wt");
	FILE *n=fopen("n.dat","wt");
	FILE *cycle=fopen("cycle.dat","wt");
	FILE *show_difference=fopen("show_difference.dat","wt");

	int i;

	float iter=50;

	float d1=1.0;
	float Error;

	float r=0.5,x,y,ceta;
	float PI=3.1415926535;
	float Umn, Vmn, Unm, Vnm, deltaT;
	float uma, vma, una, vna;
	float umb, vmb, unb, vnb;

	float xma, yma, xna, yna;

	float xmb, ymb, xmd[N], ymd[N];
	float xnb, ynb, xnd[N], ynd[N];
	
	float timestep = 3.0, step;

	step = 3;

	xma=0.0;yma=0.5;										//initialize two vortex points
	xna=0.0;yna=-0.5;										//initialize two vortex points

	deltaT=timestep/iter;									//calculat time step of iteration
	
	xmd[0]=xma; ymd[0]=yma; xnd[0]=xna; ynd[0]=yna;

	fprintf(m,"%0.6f %0.6f\n",xmd[0], ymd[0]);
	fprintf(n,"%0.6f %0.6f\n",xnd[0], ynd[0]);

	fprintf(show_difference,"%0.6f %0.6f\n",xma, yma);
	for(ceta=0;ceta<=2*PI;ceta+=2*PI/30)					//initialize cycle
	{
		x=r*cos(ceta);
		y=r*sin(ceta);
		fprintf(cycle,"%0.6f %0.6f\n",x,y);
	}

	for(i=0; i<(step * iter); i++)
	{
		Umn = cacu_U(xma, yma, xna, yna);					//calculat Umn
		Vmn = cacu_V(xma, yma, xna, yna);
		uma = cacu_u(Umn, (float)GAMn);									//calculat uma
		vma = cacu_v(Vmn, (float)GAMn);
	
		xmb = cacu_xb(xma, uma, deltaT);					//calculat xmb
		ymb = cacu_yb(yma, vma, deltaT);

		if(i==0)
		{
			fprintf(show_difference,"%0.6f %0.6f\n",xmb,ymb);
		}

		Unm = cacu_U(xna, yna, xma, yma);					//calculat Unm
		Vnm = cacu_V(xna, yna, xma, yma);	
		una = cacu_u(Unm, (float)GAMm);									//calculat una
		vna = cacu_v(Vnm, (float)GAMm);
		xnb = cacu_xb(xna, una, deltaT);					//calculat xnb
		ynb = cacu_yb(yna, vna, deltaT);
		Umn=0; Vmn=0; Unm=0; Vnm=0;

		Umn = cacu_U(xmb, ymb, xnb, ynb);					//calculat Unm
		Vmn = cacu_V(xmb, ymb, xnb, ynb);	
		umb = cacu_u(Umn, (float)GAMn);									//calculat umb
		vmb = cacu_v(Vmn, (float)GAMn);
		xmd[i+1] = cacu_xd(xma, uma, umb, deltaT);			//calculat xmd
		ymd[i+1] = cacu_yd(yma, vma, vmb, deltaT);
		
		if(i==0)
		{	
			fprintf(show_difference,"%0.6f %0.6f\n",xma,yma);
			fprintf(show_difference,"%0.6f %0.6f\n",xmb,ymb);
			fprintf(show_difference,"%0.6f %0.6f\n",xmd[i+1],ymd[i+1]);
		}

		fprintf(m,"%0.6f %0.6f\n",xmd[i+1], ymd[i+1]);

		Unm = cacu_U(xnb, ynb, xmb, ymb);
		Vnm = cacu_V(xnb, ynb, xmb, ymb);
		unb = cacu_u(Unm, (float)GAMm);
		vnb = cacu_v(Vnm, (float)GAMm);
		xnd[i+1] = cacu_xd(xna, una, unb, deltaT);
		ynd[i+1] = cacu_yd(yna, vna, vnb, deltaT);

		fprintf(n,"%0.6f %0.6f\n",xnd[i+1], ynd[i+1]);

		Error=cacu_error(xma, yma, xna, yna, d1);

		xma=xmd[i+1]; yma=ymd[i+1]; xna=xnd[i+1]; yna=ynd[i+1];
	}
	for(i=0; i<(step * iter); i++)
	{
		if(i==(iter-1)||i==(2*iter-1)||i==(3*iter-1)||i==(4*iter-1)||i==(5*iter-1)||i==(6*iter-1))
			printf("Xmd = %0.6f , Ymd = %0.6f\n",xmd[i+1], ymd[i+1]);
	}
	printf("\n");
	for(i=0; i<(step * iter); i++)
	{
		if(i==(iter-1)||i==(2*iter-1)||i==(3*iter-1)||i==(4*iter-1)||i==(5*iter-1)||i==(6*iter-1))
			printf("Xnd = %0.6f , Ynd = %0.6f\n",xnd[i+1], ynd[i+1]);
	}
}

float cacu_U(float Xm, float Ym, float Xn, float Yn)
{
	float PI=3.1415926535;
	float temp_U, tempXm, tempYm, tempXn, tempYn;
	tempXm=Xm; tempYm=Ym; tempXn=Xn; tempYn=Yn;
	temp_U=(1/(2*PI))*( (tempYm-tempYn) / ( (float)pow((tempXm-tempXn),2) + (float)pow((tempYm-tempYn),2) ));
	return temp_U;
	temp_U=0; tempXm=0; tempYm=0; tempXn=0; tempYn=0;
}
float cacu_V(float Xm, float Ym, float Xn, float Yn)
{
	float PI=3.1415926535;
	float temp_V, tempXm, tempYm, tempXn, tempYn;
	tempXm=Xm; tempYm=Ym; tempXn=Xn; tempYn=Yn;
	temp_V=(-1/(2*PI))*((tempXm-tempXn) / ( (float)pow((tempXm-tempXn),2) + (float)pow((tempYm-tempYn),2) ));
	return temp_V;
	temp_V=0; tempXm=0; tempYm=0; tempXn=0; tempYn=0;
}
float cacu_u(float U, float GAM)
{
	float temp_gam;
	float PI=3.1415926535;
	float temp_u, tempU;
	tempU=U; temp_gam=GAM;
	temp_u=temp_gam * tempU;
	return temp_u;
	temp_u=0; tempU=0; temp_gam=0;
}
float cacu_v(float V, float GAM)
{
	float temp_gam;
	float PI=3.1415926535;
	float temp_v, tempV;
	tempV=V; temp_gam=GAM;
	temp_v=temp_gam * tempV;
	return temp_v;
	temp_v=0; tempV=0; temp_gam=0;
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
float cacu_error(float X1, float Y1, float X2, float Y2, float D1)
{
	float temp_error, tempX1, tempY1, tempX2, tempY2, tempD1, tempDn;
	tempX1=X1; tempY1=Y1; tempX2=X2; tempY2=Y2; tempD1=D1;
	tempDn=(float)sqrt(( (float)pow((tempX1-tempX2),2) + (float)pow((tempY1-tempY2),2) ));
	temp_error=(tempDn-tempD1)/tempDn;
	return temp_error;
	temp_error=0; tempX1=0; tempY1=0; tempX2=0; tempY2=0; tempD1=0; tempDn=0;
}