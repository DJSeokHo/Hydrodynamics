/*********************************************************
*****    Numerical Analysis for 2D Lifting Bodies    *****
*****               2Dwing.c + 2Dwing.h              *****
**********************************************************
*****     Descriptions of the defined parameters     *****
*****                                                *****  
*****      mode: 0: Elliptic Foil                    *****  
*****            1: NACA 4digit Foil                 *****  
*****			 2: NACA Cambered Foil                     *****
*****      NP: Total panel number (should be even)   *****  
*****      Chord: Chord Length (1)                   *****  
*****      M: Maximum ordinate of mean line (Camber) *****  
*****      P: Chordwise position of maximum ordinate *****
*****      thick: Maximum Thickness                  *****
*****                                                *****  
*****   ex) NACA4415: M=0.04, P=0.4, thick=0.15      *****  
*****                                                *****  
*****      alpha: Angle of attack(degree)            ***** 
**********************************************************/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "2Dwing.h"
#define pi 3.141592

#define mode 1     
#define NP 100      
#define Chord 1.0 
#define M 0.04      
#define P 0.4
#define thick 0.12  
#define alpha 0.0  

void mkpanel(void);
void setmat(void);
void getpot(int,int,double *,double *);
void force(void);

double yt(double,double);
double yc(double, double, double);

double *Xnodp,*Ynodp,*Xcolp,*Ycolp, *Plnth, *Cp, *Vs;
double *Tx,*Ty,*Nx,*Ny;

FILE *fp_geo, *fp_pot, *fp_vel, *fp_pre ;

void main(){
  fp_geo = fopen("geo.dat","wt");
  fp_pot = fopen("pot.dat","wt");
  fp_vel = fopen("vel.dat","wt");
  fp_pre = fopen("pre.dat","wt");
  
  Xnodp=dvector(1,NP+1);  /*Node points*/
  Ynodp=dvector(1,NP+1);
  Xcolp=dvector(1,NP);    /*collocation points*/
  Ycolp=dvector(1,NP);
  Plnth=dvector(0,NP+1);  /*length of the panels*/
  Nx=dvector(1,NP);       /*normal vectors*/
  Ny=dvector(1,NP);
  Tx=dvector(1,NP);       /*tangential vectors*/
  Ty=dvector(1,NP);
  Cp=dvector(1,NP);
  Vs=dvector(1,NP);


  mkpanel();
  setmat();
  force();
  
  free_dvector(Xnodp,1,NP+1);
  free_dvector(Ynodp,1,NP+1);
  free_dvector(Xcolp,1,NP);
  free_dvector(Ycolp,1,NP);
  free_dvector(Plnth,0,NP+1);
  free_dvector(Nx,1,NP);
  free_dvector(Ny,1,NP);
  free_dvector(Tx,1,NP);
  free_dvector(Ty,1,NP);
  free_dvector(Cp,1,NP);
  free_dvector(Vs,1,NP);

  fclose(fp_geo);
  fclose(fp_pot);
  fclose(fp_vel);
  fclose(fp_pre);
}

/*******************************************************************/
/*****          Generate panels (cosine spacing)               *****/
/*******************************************************************/
void mkpanel(){
  int i;
  double nceta,arc,*xn,*yn;
  fp_geo = fopen("geo.dat","wt");

  xn=dvector(1,NP+1);
  yn=dvector(1,NP+1);
  arc=2.0*pi/NP;
  switch(mode){
/************************ Elliptic Foil ****************************/
    case 0:
    for(i=1;i<=NP+1;i++){
      nceta=(i-1)*arc;
      xn[i]=0.5*(1.+cos(nceta));
      yn[i]=0.5*thick*sin(nceta);
      Xnodp[i]=xn[i];
      Ynodp[i]=yn[i];
    }
    break;
        
/**************** NACA 4digit Symetric Foil ************************/
    case 1:
    for(i=1;i<=NP;i++){
      Xnodp[i]=Ynodp[i]=0.0;
    }
    
    for(i=1;i<=NP+1;i++){
      nceta=(i-1)*arc;
      xn[i]=0.5*(1.0+ cos(nceta));
    }
    Xnodp[1]=1.0;
    Ynodp[1]=0.0;
    yn[1]=0.0;
    yn[NP/2+1]=0.0;
    
    for(i=2;i<=NP/2;i++){
      yn[i]=yt(xn[i],thick);
      Xnodp[i]=xn[i];
      Ynodp[i]=yn[i];
    }
    for(i=NP/2+1;i<=NP+1;i++){
      Xnodp[i]=xn[i];
      Ynodp[i]=-yn[NP+2-i];
    }
    break;

/***************** NACA 4digit Cambered Foil ***********************/
    case 2:
    for(i=1;i<=NP;i++){
      Xnodp[i]=Ynodp[i]=0.0;
    }
    
    for(i=1;i<=NP+1;i++){
      nceta=(i-1)*arc;
      xn[i]=0.5*(1.0+ cos(nceta));
    }
    Xnodp[1]=1.0;
    Ynodp[1]=0.0;
    yn[1]=0.0;
    yn[NP/2+1]=0.0;
    
    for(i=2;i<=NP/2;i++){
      yn[i]=yc(xn[i],M,P)+yt(xn[i],thick);
      Xnodp[i]=xn[i];

      Ynodp[i]=yn[i];

    }
    for(i=NP/2+1;i<=NP+1;i++){
      yn[i]=yc(xn[i],M,P)-yt(xn[i],thick);
      Xnodp[i]=xn[i];
      Ynodp[i]=yn[i];

    }
    break;
/*******************************************************************/
  }

  for(i=1;i<=NP;i++){
    Xcolp[i]=(Xnodp[i+1]+Xnodp[i])*0.5;
    Ycolp[i]=(Ynodp[i+1]+Ynodp[i])*0.5;
    Plnth[i]=sqrt(pow(Xnodp[i+1]-Xnodp[i],2)+pow(Ynodp[i+1]-Ynodp[i],2));
    Nx[i]=(Ynodp[i+1]-Ynodp[i])/Plnth[i];
    Ny[i]=(Xnodp[i]-Xnodp[i+1])/Plnth[i];
    Tx[i]=-Ny[i]; /*(Xnodp[i+1]-Xnodp[i])/Plnth[i];*/
    Ty[i]= Nx[i]; /*(Ynodp[i+1]-Ynodp[i])/Plnth[i];*/
  }
  Plnth[0]=Plnth[1];
  Plnth[NP+1]=Plnth[NP];

/*******************************************************************/

  fprintf(fp_geo,"title =\"thick=%3.2f  Np=%d \"\n",thick,NP);
  fprintf(fp_geo,"zone t=\"Node point\"\n");
  for(i=1;i<=NP+1;i++)
    fprintf(fp_geo,"%3.6f %3.6f \n",Xnodp[i],Ynodp[i]);
  fprintf(fp_geo,"zone t=\"Control point\"\n");
  for(i=1; i<=NP ;i++)
  fprintf(fp_geo,"%3.6f %3.6f \n",Xcolp[i],Ycolp[i]);

  fclose(fp_geo);
}

/**********************************/
/***** Solve the Algebraic Eq *****/
/**********************************/
void setmat(void){
  int i,j,self=0,*indx;
  double Uxinf,Uyinf,dx1,dx2,d,rx,ry,kceta,KUTTA,Cl;
  double *phi,*RHS,*Source,*Dblet,*Sigma,*indvs;
  double **Aij,**Bij,**cd;

  phi=dvector(0,NP+1);    /*induced velocity potential*/
  RHS=dvector(1,NP);
  Source=dvector(1,NP);
  Dblet=dvector(1,NP);
  cd=dmatrix(1,NP,1,NP);  /*croneca delta*/
  Aij=dmatrix(1,NP,1,NP);
  Bij=dmatrix(1,NP,1,NP);
  Sigma=dvector(1,NP);    /*source strength*/
  indvs=dvector(1,NP);
  indx=ivector(1,NP);
  
  Uxinf=cos(alpha*pi/180.0); /*free stream velocity*/
  Uyinf=sin(alpha*pi/180.0);
  for(i=1;i<=NP;i++){
    RHS[i]=0.0;
    Sigma[i]=(Nx[i]*Uxinf+Ny[i]*Uyinf); /*source strength*/
  }
  for(i=1;i<=NP;i++){
    self=i;
    getpot(NP,self,Dblet,Source);
    for(j=1;j<=NP;j++){
      cd[j][i]=0.0;
      cd[j][j]=1.0;
      Aij[j][i]=cd[j][i]-Dblet[j];
      Bij[j][i]=Source[j];
      RHS[j]=RHS[j]-Sigma[i]*Bij[j][i];
    }
  }
  /**KUTTA Condition**/
  rx=Xcolp[1]-Xcolp[NP];
  ry=Ycolp[1]-Ycolp[NP];
  for(j=1;j<=NP;j++){
    kceta=atan2(Ycolp[j],Xcolp[j]-1.0);
    if(kceta<0.0) {KUTTA=(kceta+2.0*pi)/(2.0*pi);}
    else {KUTTA=kceta/(2.0*pi);}
    Aij[j][1]=Aij[j][1]+KUTTA;
    Aij[j][NP]=Aij[j][NP]-KUTTA;
    RHS[j]=RHS[j]-(Uxinf*rx+Uyinf*ry)*KUTTA;
  }
  
  ludcmp(Aij,NP,indx,&d);
  lubksb(Aij,NP,indx,RHS);
  
  for(i=1;i<=NP;i++) phi[i]=RHS[i];
  phi[0]=phi[1];
  phi[NP+1]=phi[NP];
  Cl=2.0*(phi[1]-phi[NP]);

  /** 2nd order Polynormial(central)**/
  for(i=1;i<=NP;i++){
    if(i>=1 && i<=3){
      dx1=0.5*(Plnth[i]+Plnth[i+1]);
      dx2=0.5*(Plnth[i+1]+Plnth[i+2]);
      indvs[i]=-(dx1*dx1*phi[i+2]-
		 (dx1+dx2)*(dx1+dx2)*phi[i+1]+
		 (2.*dx1*dx2+dx2*dx2)*phi[i])/
	(dx1*dx2*(dx1+dx2));
    }
    else if(i>=NP-2 && i<=NP){
      dx1=0.5*(Plnth[i]+Plnth[i-1]);
      dx2=0.5*(Plnth[i-1]+Plnth[i-2]);
      indvs[i]=(dx1*dx1*phi[i-2]-
		(dx1+dx2)*(dx1+dx2)*phi[i-1]+
		(2.*dx1*dx2+dx2*dx2)*phi[i])/
	(dx1*dx2*(dx1+dx2));
    }
    else{
      dx1=0.5*(Plnth[i]+Plnth[i-1]);
      dx2=0.5*(Plnth[i]+Plnth[i+1]);
      indvs[i]=(dx1*dx1*(phi[i+1]-phi[i])+
		dx2*dx2*(phi[i]-phi[i-1]))/
	(dx1*dx2*(dx1+dx2));
    }
    Vs[i]=-(Uxinf*Tx[i]+Uyinf*Ty[i]+indvs[i]);
    Vs[1]=Vs[NP]=0.0;
    Cp[i]=1.0-(Vs[i]*Vs[i]);
  }
  fprintf(fp_pot,"title=\"t=%2.2f %dPANEL alpha=%2.1f \"\n",thick,NP,alpha);
  fprintf(fp_vel,"title=\"t=%2.2f %dPANEL alpha=%2.1f \"\n",thick,NP,alpha);
  fprintf(fp_pre,"title=\"t=%2.2f %dPANEL alpha=%2.1f \"\n",thick,NP,alpha);
  fprintf(fp_pot,"variables=\"x\",\"y\"\n");
  fprintf(fp_vel,"variables=\"x\",\"y\"\n");
  fprintf(fp_pre,"variables=\"x\",\"y\"\n");
  fprintf(fp_pot,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl);
  fprintf(fp_vel,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl);
  fprintf(fp_pre,"zone t=\"h/c=inf  Cl=%3.5f\"\n",Cl);
  for(i=1;i<=NP;i++){
    fprintf(fp_pot,"%3.6f  %3.6f  \n",Xcolp[i],phi[i]);
    fprintf(fp_vel,"%3.6f  %3.6f  \n",Xcolp[i],fabs(Vs[i]));
    fprintf(fp_pre,"%3.6f  %3.6f  \n",Xcolp[i],-Cp[i]);
  }

  free_dvector(phi,0,NP+1);
  free_dvector(RHS,1,NP);
  free_dvector(Source,1,NP);
  free_dvector(Dblet,1,NP);
  free_dmatrix(cd,1,NP,1,NP);
  free_dmatrix(Aij,1,NP, 1,NP);
  free_dmatrix(Bij,1,NP, 1,NP);
  free_dvector(Sigma,1,NP);
  free_dvector(indvs,1,NP);
  free_ivector(indx,1,NP);
}


/******************************************/
/***** Computes the induced potential *****/
/******************************************/
void getpot(int NNP,int k,double *dblet,double *source){
  int j;
  double tx,ty,xn1,yn1,xn2,yn2,xc,yc,pl;
  double	l1,l2,twr1,twr2,h,theta;
  xn1=Xnodp[k];
  yn1=Ynodp[k];
  xn2=Xnodp[k+1];
  yn2=Ynodp[k+1];
  pl=Plnth[k];
  tx=(xn2-xn1)/pl;
  ty=(yn2-yn1)/pl;
  for(j=1;j<=NNP;j++){
    xc=Xcolp[j];
    yc=Ycolp[j];
    l1=(xc-xn1)*tx+(yc-yn1)*ty;
    l2=-((xc-xn2)*tx+(yc-yn2)*ty);
    h=-((xc-xn1)*ty-(yc-yn1)*tx);
    twr1=(xc-xn1)*(xc-xn1)+(yc-yn1)*(yc-yn1); /*squre(r1)*/
    twr2=(xc-xn2)*(xc-xn2)+(yc-yn2)*(yc-yn2); /*squre(r2)*/
    theta=atan2(fabs(h)*pl,twr1-pl*l1);
    if (j==k){
      dblet[j]=0.5;
      source[j]=0.0;
    }
    else
      dblet[j]=-h/fabs(h)*theta/(2.0*pi);
      source[j]=0.5*((pl-l1)*log(twr2)+l1*log(twr1)-2.0*pl+2.0*fabs(h)*theta)/(2.0*pi);
  }
}



/*******************************************************************/
/******       Compute forces (CL, CD)                           *****/
/*******************************************************************/
void force(void){
  int i;
  double Fy,Fx,Fxd,Fd,cf;
  double CL,CD,CLf,CDf;

  Fy=Fx=Fd=Fxd=0.0;
  CL=CLf=CD=CDf=0.0;
  cf=0.0035;

  for(i=1; i<=NP; i++){
    Fy  +=   -Cp[i]*Ny[i]*Plnth[i];
    Fx  +=   -Cp[i]*Nx[i]*Plnth[i];
    Fd  += cf*Vs[i]*Vs[i]*Plnth[i];
  }
    Fxd=Fx+Fd;
  
  CL = Fy*cos(alpha*pi/180.0) - Fx*sin(alpha*pi/180.0);
  CD = Fx*cos(alpha*pi/180.0) + Fy*sin(alpha*pi/180.0);
 
  CLf = Fy*cos(alpha*pi/180.0) - Fxd*sin(alpha*pi/180.0);
  CDf = Fxd*cos(alpha*pi/180.0) + Fy*sin(alpha*pi/180.0);

printf("\n==========================================================\n");
printf(" Ang(deg)     CL         CD          CLf         CDf\n");
printf("----------------------------------------------------------\n");
printf(" %3.1f       %3.5f    %3.5f    %3.5f      %3.5f", alpha,CL,CD,CLf,CDf);
printf("\n----------------------------------------------------------\n");

}


/*******************************************************************/
/*****            NACA foil (Thickness Distribution)           *****/
/*******************************************************************/
double yt(double x, double t){
  double y=0.0;
  if(x==1.0) return y=0.0;
  if(x==0.0) return y=0.0;
  else {
    y=(t/0.2)*(0.2969*pow(x,0.5)-0.126*x-0.3537*x*x
	       +0.2843*pow(x,3.)-0.1015*pow(x,4.));
    return y;
  }
}
/*******************************************************************/
/*****              4digit NACA foil (Meanline)                *****/
/*******************************************************************/
double yc(double x, double m, double p)
{
  double y;
  if(x<=p) return y=m*(2*p*x-x*x)/(p*p);
  else return y=m*((1-2*p)+2*p*x-x*x)/(1-2*p+p*p);
}




