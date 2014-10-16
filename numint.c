#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jgm3.h"
#include "vecmat.h"
#include "refsys.h"

double mjd0;

void rkck(double y[], double dydx[], int n, double x, double h, double yout[],double yerr[], void (*derivs)(double, double [], double []))
{
  int i;
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
  
  ak2=(double *) malloc(sizeof(double)*n);
  ak3=(double *) malloc(sizeof(double)*n);
  ak4=(double *) malloc(sizeof(double)*n);
  ak5=(double *) malloc(sizeof(double)*n);
  ak6=(double *) malloc(sizeof(double)*n);
  ytemp=(double *) malloc(sizeof(double)*n);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  (*derivs)(x+a2*h,ytemp,ak2);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*derivs)(x+a3*h,ytemp,ak3);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*derivs)(x+a4*h,ytemp,ak4);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*derivs)(x+a5*h,ytemp,ak5);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (*derivs)(x+a6*h,ytemp,ak6);
  for (i=0;i<n;i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=0;i<n;i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);

  free(ytemp);
  free(ak6);
  free(ak5);
  free(ak4);
  free(ak3);
  free(ak2);

  return;
}

void accel(double t,double r[],double drdt[])
{
  int i;
  double mjd;
  double e[3][3],rin[3],a[3];

  // Get time
  mjd=mjd0+t/86400.0;

  // Compute transformation
  icrs_to_itrs(mjd,e);

  // Input vector
  for (i=0;i<3;i++)
    rin[i]=r[i];

  // Compute geopotential
  jgm3_geopotential(rin,e,20,20,a);

  // Derivatives
  for (i=0;i<3;i++) {
    drdt[i]=r[i+3];
    drdt[i+3]=a[i];
  }

  return;
}

int main(int argc,char *argv[])
{
  int i,j,n=6;
  double r[6],drdt[6],rout[6],rerr[6];
  double t,dt;
  
  mjd0=56946.75000;
  r[0]=5315.36;
  r[1]=4225.90;
  r[2]=-220.61;
  r[3]=-3.109;
  r[4]=3.604;
  r[5]=-6.006;
  t=0.0;
  dt=10.0;
  
  for (i=0;i<1000;i++) {
    accel(t,r,drdt);
    rkck(r,drdt,n,t,dt,rout,rerr,accel);
    
    for (j=0;j<n;j++)
      r[j]=rout[j];
    printf("%14.8lf %f %f %f %f %f %f\n",mjd0+t/86400.0,r[0],r[1],r[2],r[3],r[4],r[5]);
    t+=dt;
  }
  

  return 0;
}
