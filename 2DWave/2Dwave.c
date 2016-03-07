#include <stdio.h>
#include <math.h>

#define U 1
#define g 9.8
#define PI 3.1415926535

void wavedrag();

main()
{
	wavedrag();
	return 0;
}

void wavedrag()
{
	FILE *wave=fopen("wave.dat","wt");
	FILE *length=fopen("length.dat","wt");
	double D, F;
	double l;
	double wavelength;

	int i;
	wavelength = (2*PI*pow(U,2)) / g;
	printf("wavelength = %0.6lf\n\n",wavelength);
	
	printf("When l approximate (n * wavelength),the wave drag approximate zero:\n");  
	for(i=1;i<=1.6*1000;i++)
	{
		F=(double)i/1000;
		l=pow(U,2) / (pow(F,2) * g);
		D = pow(sin(l*g*0.5),2);
		fprintf(wave,"%0.6lf  %0.6lf\n",F,D);
		fprintf(length,"%0.6lf  %0.6lf  %0.6lf\n",F,D,l);
		if(D<0.0001)
		{
			printf( "l approximate %0.6lf   D approximate %0.6lf\n", l, D );
		}
	}

	printf("\nWhen l approximate ((n+0.5) * wavelength),the wave drag approximate one:\n");  
	for(i=1;i<=1.6*1000;i++)
	{
		F=(double)i/1000;
		l=pow(U,2) / (pow(F,2) * g);
		D = pow(sin(l*g*0.5),2);
		if(D>0.99950 && l>0.32)
		{
			printf( "l approximate %0.6lf   D approximate %0.6lf\n", l, D );
		}
	}
	return 0;
}


