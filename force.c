#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jgm3.h"
#include "force.h"
#include "vecmat.h"
#include "constants.h"
#include "refsys.h"

double density(double r[3],double s[3])
{
  int i,n_prm=3,n=50,ih;
  double ra,de,ra_lag=0.523599; // Lag in RA [rad]
  double h0,hmin=100.0,hmax=1000.0; // Model limits;
  double u[3],cospsi2,hsmin,hsmax,rho,rhomin,rhomax;
  double h[]={100.0,120.0,130.0,140.0,150.0,160.0,170.0,180.0,190.0,200.0,
	      210.0,220.0,230.0,240.0,250.0,260.0,270.0,280.0,290.0,300.0,
	      320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,
	      520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,
	      720.0,740.0,760.0,780.0,800.0,840.0,880.0,920.0,960.0,1000.0};
  double rhohmin[]={4.974e+02,2.490e+01,8.377e+00,3.899e+00,2.122e+00,1.263e+00,
		    8.008e-01,5.283e-01,3.617e-01,2.557e-01,1.839e-01,1.341e-01,
		    9.949e-02,7.488e-02,5.709e-02,4.403e-02,3.430e-02,2.697e-02,
		    2.139e-02,1.708e-02,1.099e-02,7.214e-03,4.824e-03,3.274e-03,
		    2.249e-03,1.558e-03,1.091e-03,7.701e-04,5.474e-04,3.916e-04,
		    2.819e-04,2.042e-04,1.488e-04,1.092e-04,8.070e-05,6.012e-05,
		    4.519e-05,3.430e-05,2.632e-05,2.043e-05,1.607e-05,1.281e-05,
		    1.036e-05,8.496e-06,7.069e-06,4.680e-06,3.200e-06,2.210e-06,
		    1.560e-06,1.150e-06};
  double rhohmax[]={4.974e+02,2.490e+01,8.710e+00,4.059e+00,2.215e+00,1.344e+00,
		    8.758e-01,6.010e-01,4.297e-01,3.162e-01,2.396e-01,1.853e-01,
		    1.455e-01,1.157e-01,9.308e-02,7.555e-02,6.182e-02,5.095e-02,
		    4.226e-02,3.526e-02,2.511e-02,1.819e-02,1.337e-02,9.955e-03,
		    7.492e-03,5.684e-03,4.355e-03,3.362e-03,2.612e-03,2.042e-03,
		    1.605e-03,1.267e-03,1.005e-03,7.997e-04,6.390e-04,5.123e-04,
		    4.121e-04,3.325e-04,2.691e-04,2.185e-04,1.779e-04,1.452e-04,
		    1.190e-04,9.776e-05,8.059e-05,5.741e-05,4.210e-05,3.130e-05,
		    2.360e-05,1.810e-05};

  // Geodetic height of satellite
  h0=geodetic_height(r);

  // Compute density
  if (h0>=hmin && h0<=hmax) {
    // Compute solar subpoint
    ra=atan2(s[1],s[0]);
    de=atan2(s[2],sqrt(s[0]*s[0]+s[1]*s[1]));

    // Unit vector towards lagged apex
    u[0]=cos(de)*cos(ra+ra_lag);
    u[1]=cos(de)*sin(ra+ra_lag);
    u[2]=sin(de);

    // Angle between satellite and apex
    cospsi2=0.5+0.5*dot_product(r,u)/magnitude(r);

    // Locate segment
    ih=0;
    for (i=0;i<n-1;i++) {
      if (h0>=h[i] && h0<h[i+1]) {
	ih=i;
	break;
      }
    }

    // Compute scale heights
    hsmin=(h[ih]-h[ih+1])/log(rhohmin[ih+1]/rhohmin[ih]);
    hsmax=(h[ih]-h[ih+1])/log(rhohmax[ih+1]/rhohmax[ih]);
    rhomin=rhohmin[ih]*exp((h[ih]-h0)/hsmin);
    rhomax=rhohmax[ih]*exp((h[ih]-h0)/hsmax);

    // Density
    rho=rhomin+(rhomax-rhomin)*pow(cospsi2,n_prm);
  } else {
    rho=0.0;
  }
		    
  return rho;
}


// Drag
void drag(double r[3],double v[3],double s[3],double e[3][3],double area,double mass,double cd,double a[3])
{
  int i;
  double rtod[3],vtod[3],et[3][3],vrel[3],rv,vc[3],atod[3],fac;
  double omega[3]={0.0,0.0,7.292114992e-5},rho;

  // Convert to EME2000 FIXME
  vector_multiply(e,r,rtod);
  vector_multiply(e,v,vtod);

  // Relative velocity
  cross_product(omega,rtod,vc);
  for (i=0;i<3;i++) 
    vrel[i]=vtod[i]-vc[i];
  rv=magnitude(vrel);

  // Density
  rho=density(r,s);

  // Acceleration
  fac=-0.5*cd*(area/mass)*rho*rv;
  for (i=0;i<3;i++)
    atod[i]=fac*vrel[i];

  // Back to ICRS
  matrix_transpose(e,et);
  vector_multiply(et,atod,a);

  return;
}

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
  rs=dot_product(r,e);

  // 
  for (i=0;i<3;i++)
    d[i]=r[i]-e[i]*rs;
  rd=magnitude(d);

  // Illumination fraction
  return ((rs>0 || rd>XKMPER) ? 1.0 : 0.0);
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

