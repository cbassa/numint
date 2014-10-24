#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "vecmat.h"
#include "constants.h"
#include "series.h"
#include "jpleph.h"

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}

// Greenwich Mean Sidereal Time
double gmst(double mjd)
{
  double t,gmst;

  t=(mjd-51544.5)/36525.0;

  gmst=modulo(280.46061837+360.98564736629*(mjd-51544.5)+t*t*(0.000387933-t/38710000),360.0)*D2R;

  return gmst;
}

// Nutation
void nutation(double mjd,double *dpsi,double *deps,double *eps)
{
  int i;
  double t,t2,t3;
  double d,m,m1,f,n;
  double arg;

  // Julian centuries
  t=(mjd-51544.5)/36525.0;
  t2=t*t;
  t3=t2*t;

  // Angles
  d=modulo(297.85036+445267.111480*t-0.0019142*t2+t3/189474.0,360.0)*D2R;
  m=modulo(357.52772+35999.050340*t-0.0001603*t2-t3/300000.0,360.0)*D2R;
  m1=modulo(134.96298+477198.867398*t+0.0086972*t2+t3/56250.0,360.0)*D2R;
  f=modulo(93.27191+483202.017538*t-0.0036825*t2+t3/327270.0,360.0)*D2R;
  n=modulo(125.04452-1934.136261*t+0.0020708*t2+t3/450000.0,360.0)*D2R;
  
  // Compute sums
  *dpsi=0.0;
  *deps=0.0;
  for (i=0;i<106;i++) {
    arg=nut[i].nd*d+nut[i].nm*m+nut[i].nm1*m1+nut[i].nf*f+nut[i].nn*n;
    *dpsi+=(nut[i].s+nut[i].st*t)*sin(arg);
    *deps+=(nut[i].c+nut[i].ct*t)*cos(arg);
  }
  *dpsi*=0.0001/3600*D2R;
  *deps*=0.0001/3600*D2R;
  *eps=-46.8150*t-0.00059*t2+0.001813*t3;
  *eps=(23.4392911+ *eps/3600.0)*D2R;

  return;
}

// Precession
void precess(double mjd0,double mjd,double *zeta,double *z,double *theta)
{
  double t0,t;

  // Time in centuries
  t0=(mjd0-51544.5)/36525.0;
  t=(mjd-mjd0)/36525.0;

  // Precession angles
  *zeta=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  *zeta+=(0.30188-0.000344*t0)*t*t+0.017998*t*t*t;
  *zeta*=D2R/3600.0;
  *z=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  *z+=(1.09468+0.000066*t0)*t*t+0.018203*t*t*t;
  *z*=D2R/3600.0;
  *theta=(2004.3109-0.85330*t0-0.000217*t0*t0)*t;
  *theta+=-(0.42665+0.000217*t0)*t*t-0.041833*t*t*t;
  *theta*=D2R/3600.0;
  
  return;
}

// ECI to ECEF
void eci_to_ecef(double mjd,double g[3][3])
{
  int i,j;
  double h,dpsi,deps,eps;

  // Nutation
  nutation(mjd,&dpsi,&deps,&eps);

  // Earth orientation
  h=gmst(mjd);
  identity_matrix(g);
  rotate_z(h+dpsi*cos(eps+deps),g);

  return;
}

// ICRS to EME conversion
void icrs_to_eme(double mjd,double a[3][3])
{
  int i,j;
  double dpsi,deps,eps,z,theta,zeta,h;
  double p[3][3],n[3][3];

  // Precession
  precess(51544.5,mjd,&zeta,&z,&theta);
  identity_matrix(p);
  rotate_z(-zeta,p);
  rotate_y(theta,p);
  rotate_z(-z,p);

  // Nutation
  nutation(mjd,&dpsi,&deps,&eps);
  identity_matrix(n);
  rotate_x(eps,n);
  rotate_z(-dpsi,n);
  rotate_x(-eps-deps,n);

  // Multiply matrices (left to right)
  matrix_multiply(n,p,a);

  return;
}

// ICRS to ITRS conversion
void icrs_to_itrs(double mjd,double a[3][3])
{
  int i,j;
  double dpsi,deps,eps,z,theta,zeta,h;
  double p[3][3],n[3][3],g[3][3],c[3][3],b[3][3];

  // Precession
  precess(51544.5,mjd,&zeta,&z,&theta);
  identity_matrix(p);
  rotate_z(-zeta,p);
  rotate_y(theta,p);
  rotate_z(-z,p);

  // Nutation
  nutation(mjd,&dpsi,&deps,&eps);
  identity_matrix(n);
  rotate_x(eps,n);
  rotate_z(-dpsi,n);
  rotate_x(-eps-deps,n);

  // Earth orientation
  h=gmst(mjd);
  identity_matrix(g);
  rotate_z(h+dpsi*cos(eps+deps),g);

  // Polar motion
  identity_matrix(c);

  // Multiply matrices (left to right)
  matrix_multiply(c,g,a);
  matrix_multiply(a,n,b);
  matrix_multiply(b,p,a);

  return;
}

// Simple solar position
void sun_position(double mjd,void *ephem,double r_sun[3])
{
  int err_code;
  int i;
  double rv[6];

  err_code=jpl_pleph(ephem,mjd+2400000.5,11,3,rv,0);

  for (i=0;i<3;i++)
    r_sun[i]=rv[i]*XKMPAU;

  return;
}

// Moon position
void moon_position(double mjd,void *ephem,double r_moon[3])
{
  int err_code;
  int i;
  double rv[6];

  err_code=jpl_pleph(ephem,mjd+2400000.5,10,3,rv,0);

  for (i=0;i<3;i++)
    r_moon[i]=rv[i]*XKMPAU;

  return;
}

// Geodetic height
double geodetic_height(double r[3])
{
  double h;

  // TO DO: COMPUTE ACTUAL GEODETIC HEIGHT

  h=magnitude(r)-XKMPER;

  return h;
}
