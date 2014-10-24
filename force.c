#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jgm3.h"
#include "force.h"
#include "vecmat.h"
#include "constants.h"

// Solar radiation pressure
// a:  area [m^2]
// m:  mass [kg]
// cr: Radiation pressure coefficient
void srp(double r[3],double s[3],double area,double mass,double cr,double a[3])
{
  int i;
  double d[3],rd,fac;
  
  // Offset from satellite to Sun
  for (i=0;i<3;i++) 
    d[i]=r[i]-s[i];

  // Distance
  rd=magnitude(d);

  // Factor
  fac=cr*area/mass*P0*(XKMPAU*XKMPAU)*pow(rd,-3);

  // Accelerations
  for (i=0;i<3;i++)
    a[i]=fac*d[i];
  
  return;
}

// Illumination
double illumination(double r[3],double s[3])
{
  int i;
  double e[3],rs,d[3],rd,a;

  // Distance
  rs=magnitude(s);

  // Unit vector
  for (i=0;i<3;i++)
    e[i]=s[i]/rs;

  // Projection
  rs=dot_product(r,s);

  // 
  for (i=0;i<3;i++)
    d[i]=r[i]-e[i]*rs;
  rd=magnitude(d);

  // Illumination fraction
  return ((d>0 || rd>XKMPER) ? 1.0 : 0.0);
}

// Third body accelaration
void third_body(double r[3],double s[3],double gm,double a[3])
{    
  int i;
  double d[3],rd,rs;
  
  // Offset from of satellite to third body
  for (i=0;i<3;i++)
    d[i]=r[i]-s[i];
  
  // Distances
  rd=magnitude(d);
  rs=magnitude(s);

  // Acceleration 
  for (i=0;i<3;i++)
    a[i]=-gm*(d[i]*pow(rd,-3)+s[i]*pow(rs,-3));

  return;
}


// JGM-3 geopotential
void jgm3_geopotential(double r[3],double e[3][3],int nmax,int mmax,double a[3])
{
  int i,n,m;
  double r_bf[3],rho,r_sqr,r0[3],a_bf[3],et[3][3];
  double v[nmax+2][nmax+2],w[nmax+2][nmax+2];
  double c,s,fac;

  // Apply transformation to earth fixed
  vector_multiply(e,r,r_bf);

  // Radii
  r_sqr=dot_product(r_bf,r_bf);
  rho=R_JGM3*R_JGM3/r_sqr;

  // Normalize coordinates
  for (i=0;i<3;i++)
    r0[i]=r_bf[i]*R_JGM3/r_sqr;

  // Compute zonal terms
  v[0][0]=R_JGM3/magnitude(r_bf);
  w[0][0]=0.0;
  v[1][0]=r0[2]*v[0][0];
  w[1][0]=0.0;
  for (n=2;n<=nmax+1;n++) {
    v[n][0]=((2*n-1)*r0[2]*v[n-1][0]-(n-1)*rho*v[n-2][0])/n;
    w[n][0]=0.0;
  }

  // Compute tesseral and sectorial terms
  for (m=1;m<=mmax+1;m++) {
    v[m][m]=(2*m-1)*(r0[0]*v[m-1][m-1]-r0[1]*w[m-1][m-1]);
    w[m][m]=(2*m-1)*(r0[0]*w[m-1][m-1]+r0[1]*v[m-1][m-1]);

    if (m<=nmax) {
      v[m+1][m]=(2*m+1)*r0[2]*v[m][m];
      w[m+1][m]=(2*m+1)*r0[2]*w[m][m];
    }

    for (n=m+2;n<=nmax+1;n++) {
      v[n][m]=((2*n-1)*r0[2]*v[n-1][m]-(n+m-1)*rho*v[n-2][m])/(n-m);
      w[n][m]=((2*n-1)*r0[2]*w[n-1][m]-(n+m-1)*rho*w[n-2][m])/(n-m);
    }
  }

  // Compute accelerations
  a_bf[0]=a_bf[1]=a_bf[2]=0.0;
  for (m=0;m<=mmax;m++) {
    for (n=m;n<=nmax;n++) {
      if (m==0) {
	c=CS_JGM3[n][0];
	a_bf[0]-=c*v[n+1][1];
	a_bf[1]-=c*w[n+1][1];
	a_bf[2]-=(n+1)*c*v[n+1][0];
      } else {
	c=CS_JGM3[n][m];
	s=CS_JGM3[m-1][n];
	fac=0.5*(n-m+1)*(n-m+2);
	a_bf[0]+=0.5*(-c*v[n+1][m+1]-s*w[n+1][m+1])+fac*(+c*v[n+1][m-1]+s*w[n+1][m-1]);
	a_bf[1]+=0.5*(-c*w[n+1][m+1]+s*v[n+1][m+1])+fac*(-c*w[n+1][m-1]+s*v[n+1][m-1]);
	a_bf[2]+=(n-m+1)*(-c*v[n+1][m]-s*w[n+1][m]);
      }
    }
  }

  // Multicative factor to km s^-2
  fac=GM_JGM3/(R_JGM3*R_JGM3);
  for (i=0;i<3;i++)
    a_bf[i]*=fac;

  // ITRS to ICRS transformation
  matrix_transpose(e,et);
  vector_multiply(et,a_bf,a);

  return;
}

