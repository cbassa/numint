#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "vecmat.h"
#include "constants.h"
#include "series.h"

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
void sun_position(double mjd,double r_sun[3])
{
  double jd,t,l0,m,e,c,r;
  double n,s,ecl,ra,de;

  jd=mjd+2400000.5;
  t=(jd-2451545.0)/36525.0;
  l0=modulo(280.46646+t*(36000.76983+t*0.0003032),360.0)*D2R;
  m=modulo(357.52911+t*(35999.05029-t*0.0001537),360.0)*D2R;
  e=0.016708634+t*(-0.000042037-t*0.0000001267);
  c=(1.914602+t*(-0.004817-t*0.000014))*sin(m)*D2R;
  c+=(0.019993-0.000101*t)*sin(2.0*m)*D2R;
  c+=0.000289*sin(3.0*m)*D2R;

  r=1.000001018*(1.0-e*e)/(1.0+e*cos(m+c));
  n=modulo(125.04-1934.136*t,360.0)*D2R;
  s=l0+c+(-0.00569-0.00478*sin(n))*D2R;
  ecl=(23.43929111+(-46.8150*t-0.00059*t*t+0.001813*t*t*t)/3600.0+0.00256*cos(n))*D2R;

  ra=atan2(cos(ecl)*sin(s),cos(s));
  de=asin(sin(ecl)*sin(s));

  r_sun[0]=r*cos(de)*cos(ra)*XKMPAU;
  r_sun[1]=r*cos(de)*sin(ra)*XKMPAU;
  r_sun[2]=r*sin(de)*XKMPAU;

  return;
}

// Moon position
void moon_position(double mjd,double r_moon[3])
{
  int i;
  double t,t2,t3,t4;
  double l1,d,m,m1,f,a1,a2,a3,e,ef;
  double suml,sumb,sumr,arglr,argb;
  double l,b,r,ra,de,eps,deps,dpsi;

  // Julian Centuries
  t=(mjd-51544.5)/36525.0;

  // Powers of t
  t2=t*t;
  t3=t2*t;
  t4=t3*t;

  // angles
  l1=modulo(218.3164477+481267.88123421*t-0.0015786*t2+t3/538841.0-t4/65194000.0,360.0)*D2R;
  d=modulo(297.8501921+445267.1114034*t-0.0018819*t2+t3/545868.0-t4/113065000.0,360.0)*D2R;
  m=modulo(357.5291092+35999.0502909*t-0.0001536*t2+t3/24490000.0,360.0)*D2R;
  m1=modulo(134.9633964+477198.8675055*t+0.0087417*t2+t3/69699.0-t4/14712000.0,360.0)*D2R;
  f=modulo(93.2720950+483202.0175233*t-0.0036539*t2-t3/3526000.0+t4/86331000.0,360.0)*D2R;
  a1=modulo(119.75+131.849*t,360.0)*D2R;
  a2=modulo(53.09+479264.290*t,360.0)*D2R;
  a3=modulo(313.45+481266.484*t,360.0)*D2R;
  e=1.0-0.002516*t-0.0000074*t2;

  // Compute sums
  for (i=0,suml=sumb=sumr=0.0;i<60;i++) {
    // Arguments
    arglr=clr[i].nd*d+clr[i].nm*m+clr[i].nm1*m1+clr[i].nf*f;
    argb=cb[i].nd*d+cb[i].nm*m+cb[i].nm1*m1+cb[i].nf*f;
    
    // E multiplication factor
    if (abs(clr[i].nm)==1)
      ef=e;
    else if (abs(clr[i].nm)==2)
      ef=e*e;
    else
      ef=1.0;

    // Sums
    suml+=clr[i].sa*sin(arglr)*ef;
    sumr+=clr[i].ca*cos(arglr)*ef;

    // E multiplication factor
    if (abs(cb[i].nm)==1)
      ef=e;
    else if (abs(cb[i].nm)==2)
      ef=e*e;
    else
      ef=1.0;

    // Sums
    sumb+=cb[i].sa*sin(argb)*ef;
  }

  // Additives
  suml+=3958*sin(a1)+1962*sin(l1-f)+318*sin(a2);
  sumb+=-2235*sin(l1)+382*sin(a3)+175*sin(a1-f)+175*sin(a1+f)+127*sin(l1-m1)-115*sin(l1+m1);

  // Ecliptic longitude, latitude and distance
  l=modulo(l1*R2D+suml/1000000.0,360.0)*D2R;
  b=sumb/1000000.0*D2R;
  r=385000.56+sumr/1000.0;

  // Equatorial
  nutation(mjd,&dpsi,&deps,&eps);
  eps+=deps;
  l+=dpsi;
  ra=modulo(atan2(sin(l)*cos(eps)-tan(b)*sin(eps),cos(l)),2.0*M_PI);
  de=asin(sin(b)*cos(eps)+cos(b)*sin(eps)*sin(l));

  // Position
  r_moon[0]=r*cos(de)*cos(ra);
  r_moon[1]=r*cos(de)*sin(ra);
  r_moon[2]=r*sin(de);

  return;
}
