#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define RnBase 10000

void ATTC(void);
void ITTC(void);
void Hughes(void);

main()
{
	ITTC();
	Hughes();
	ATTC();
	return 0;
}




void ATTC(void)
{
	FILE *attc=fopen("ATTC.dat","wt");
	int i,j,count;
	double c,z,yvalue,x;
	double ATTCRn,ATTCCf,ATTCdenominator,ATTClog;
	double temp[1000],avg;
	fprintf( attc, "Title=“ATTC”\n");
	fprintf( attc, "Variables=“X”,“Y”\n");
	fprintf( attc, "Zone I=20\n\n");
	for( j=50;j<=1000;j+=50){
		ATTCRn = (double)j*RnBase;
		count = 0;
		avg = 0;
		printf("No.%d\n",j/50);
		for(i=1;i<=100000;i++){
			x = (double)i/1000000;
			c = log10(ATTCRn*x);
			z = 0.242 / sqrt(x);
			yvalue = z - c;
			yvalue = fabs( yvalue );

			if( yvalue <= 0.0005 )
			{
				temp[j/50-1]=x;
				avg+=temp[j/50-1];
				count++;
	//			printf( "good\n");
			}
		}
		ATTCCf = (double)avg/(double)count;
		fprintf(attc,"%0.6lf  ",ATTCRn);
		fprintf(attc,"%0.9lf\n",ATTCCf);
	}
}

void ITTC(void)
{
	FILE *ittc = fopen("ITTC.dat","wt");
	int i;
	double ITTCRn,ITTCCf,ITTCdenominator;
	fprintf( ittc, "Title=“ITTC”\n");
	fprintf( ittc, "Variables=“X”,“Y”\n");
	fprintf( ittc, "Zone I=20\n\n");
	for( i=0;i<=1000;i+=50 ){
		ITTCRn = (double)i*RnBase;
		ITTCdenominator = log10(ITTCRn)-2;
		ITTCCf = 0.075 / pow(ITTCdenominator,2);
		if(i!=0){
			fprintf(ittc,"%0.6lf  ",ITTCRn);
			fprintf(ittc,"%0.6lf\n",ITTCCf);
		}
	}
}

void Hughes(void)
{
	FILE *hughes = fopen("HUGHES.dat","wt");
	int i;
	double HughesRn,HughesCf,Hughesdenominator;
	fprintf( hughes, "\nTitle=“Hughes”\n");
	fprintf( hughes, "Variables=“X”,“Y”\n");
	fprintf( hughes, "Zone I=20\n\n");
	for( i=0;i<=1000;i+=50 ){
		HughesRn = (double)i*RnBase;
		Hughesdenominator = log10(HughesRn)-2.03;
		HughesCf = 0.066 / pow(Hughesdenominator,2);
		if(i!=0){
			fprintf(hughes,"%0.6lf  ",HughesRn);
			fprintf(hughes,"%0.6lf\n",HughesCf);
		}
	}
}

