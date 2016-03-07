#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1000
#define BN 10

float cacu_randomNum(void);
float cacu_9_2_1_randomNum(float);
void sort(float a[], int);

main()
{
	float P=0.4, Bins[N]={0,}, Bin[BN][N]={0,}, top, bot;
    int i, j, k, count[BN]={0,};

	srand((unsigned)time(NULL));

	for(i=0;i<N;i++)
    {
	//	P = cacu_9_2_1_randomNum(P);
		P =	cacu_randomNum();
		Bins[i] = P;
	}
	
	sort(Bins, N);
	
	top = 0.1; bot = 0.0;

	for(i=0;i<BN;i++)
	{
		count[i]=0;k=0;
		for(j=0;j<N;j++)
		{
			if(Bins[j]>=bot && Bins[j]<=top)
			{
				Bin[i][k]=Bins[j];
				count[i]++;
				k++;
			}
		}
		top += 0.1;
		bot += 0.1;
	}

	top = 0.1; bot = 0.0;

	for(i=0;i<10;i++)
    {
		printf("Bin %2d  :  %0.1f - %0.1f   %3d\n", i+1, bot, top, count[i]);
		top += 0.1;
		bot += 0.1;
	}

    return 0;

}

float cacu_9_2_1_randomNum(float P)
{
	float temp_p;
	temp_p=P;
	temp_p = (float)pow((1.01316+temp_p),5);
	temp_p = temp_p - (int)temp_p;
	return temp_p;
	temp_p=0;
}

float cacu_randomNum(void)
{
	float temp_q;
//	srand((unsigned)time(NULL));
	temp_q = (float)(rand()%100)/100;
	return temp_q;
}

void sort(float a[], int n) 
{ 
	int i,j; 
	float k;
	for(i=0;i<n-1;i++)
	{ 
		for(j=n-1;j>i;j--) 
		{
			if(a[j]<a[j-1])
			{ 
				k=a[j]; 
				a[j]=a[j-1]; 
				a[j-1]=k; 
			} 
		}
	} 
} 
