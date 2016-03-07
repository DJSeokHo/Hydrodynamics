#include <stdio.h>
#include <math.h>

static double per[] = {0.0 ,0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4,
                       0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 1.0};

static double thk[] = {0.0, 0.187, 0.2932, 0.4132, 0.5814, 0.8,
                       0.9274, 0.9904, 0.9924, 0.9306, 0.807, 0.622,
                       0.3754, 0.2286, 0.1496, 0.1018, 0.0};

main()
{
	FILE *testc=fopen("testc.dat","wt");
	FILE *testt=fopen("testt.dat","wt");
	double cam[20];
	double thi[20];

	double x,p=0.2,m=1,y,t=0.12;
	int i;

	for(i=0;i<17;i++)
	{
		x=per[i];

		if(x<=p)
			y=m*(2*p*x-x*x)/(p*p);
		else
			y=m*((1-2*p)+2*p*x-x*x)/(1-2*p+p*p);

		printf("%0.5lf  ",x);
		printf("%0.5lf\n",y);

		fprintf(testc, "%0.5lf  ",x);
		fprintf(testc, "%0.5lf\n",y);
	}



	for(i=0;i<17;i++)
	{
		x=per[i];
		y=(t/0.2)*(0.2969*pow(x,0.5)-0.126*x-0.3537*x*x
			+0.2843*pow(x,3.)-0.1015*pow(x,4.));
	
		fprintf(testt, "%0.5lf  ",x);
		fprintf(testt, "%0.5lf\n",y);
	}

	return 0;
}
