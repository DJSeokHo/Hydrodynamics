#include <stdio.h>

#define N 3		//...number of base point

main()
{
	int i, j;
	double Li[N], x, y=0.0, xi[N], xj[N], yi[N], fz[N], fm[N];
	double speclcase=0.0;

	x = 1.5;
	yi[0] = 1.0;
	yi[1] = 4.0;
	yi[2] = 9.0;	

	for(i=0;i<N;i++)
	{
		xi[i] = (double)i + 1.0;
		fz[i] = 1.0;
		fm[i] = 1.0;
	}

	//...Calculate Lagrance result polynomials Li
	for(i=0;i<N;i++)
	{
		for(j=0;j<=N-1;j++)
		{
			if(i!=j)
			{
				fz[i] *= (x - xi[j]);
				fm[i] *= (xi[i] - xi[j]);
				Li[i] = fz[i] / fm[i];
			}
		}
	//	printf("%d  %lf\n", i+1, Li[i]);	//..check value
		y += yi[i] * Li[i];					//..y is result
		speclcase += Li[i];					//..special case
	}
//	printf("%0.2f  %0.2f\n", y, speclcase);	//..check value
	return 0;
}
