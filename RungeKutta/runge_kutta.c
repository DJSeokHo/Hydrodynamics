#include <stdio.h>
#include <math.h>

#define N 200

float differential_function(float);
float function(float, float);

main()
{
	int i, j;
	float X[N], Y[N], sub_x, sub_y, x0, x1, y0, h, length_step, sub_n_step, n, Real_y;
	float k1, k2, a, R1, R2, err;

	FILE *fun=fopen("fun.dat","wt");
	FILE *RK_fun=fopen("RK_fun.dat","wt");
	
	printf("***********************************************************************\n");
	printf("************************* Runga Kutta method **************************\n");

	x0=0.0;
	x1=5.0;
	y0=2.0;
	h=0.01;
	R1 = 0.5;
	R2 = 0.5;
	a = 1.0;
	length_step=0.5;
	err=0.0;

	n = ((x1-x0)/length_step);	

	sub_n_step = length_step/h;

	sub_x=x0;
	sub_y=y0;
	X[0]=sub_x;
	Y[0]=sub_y;

	printf("***********************************************************************\n\n");

	for(i=0; i<=(int)n; i++)
	{
		for(j=0; j<(int)sub_n_step; j++)
		{
			k1 = differential_function(sub_x);
			k2 = differential_function(sub_x+a*h);
			sub_x = sub_x+h;
			sub_y = sub_y+(R1*k1 + R2 * k2)*h;
			if(i==0 && j==0)
			{
				printf("In sub section 1(0 ~ 0.5)\n\n");
				printf("(x0 , y0)\n");
				printf("(%0.5f , %0.5f)\n", x0, y0);
				printf("(x1 , y1)\n");
				printf("(%0.5f , %0.5f)\n", sub_x, sub_y);
			}
		}
		X[i+1] = sub_x;
		Y[i+1] = sub_y;
	}


	printf("\n****************************** Result **********************************\n");

	printf("          X      |        Y        |       real Y       |        Error   \n");

	printf("------------------------------------------------------------------------\n");

	for(i=0; i<=(int)n; i++)
	{
		Real_y = function(x0 + length_step*(float) i, y0);
		fprintf(RK_fun, "%0.5f %0.5f\n", (float)i * length_step, Real_y);
		err = fabs((Real_y-Y[i])/Real_y);
		printf("%12.5f    |    %12.5f   |   %12.5f   |   %12.5f   \n", X[i], Y[i], Real_y, err);
		fprintf(fun, "%0.5f %0.5f\n", X[i], Y[i]);
	}
	
	getch();
	return 0;
}

float differential_function( float x)
{
	float a;
	a = 2 * x - 2;
	return a;
}

float function(float x, float c)
{
	float a;
	a = x * x - 2.0 * x + c;
	return a;
}



