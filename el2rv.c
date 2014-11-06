#include <stdio.h>
#include <math.h>

#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define GM_EARTH 398600.4415 // [km^3/s^2]; JGM3

int main(int argc,char *argv[])
{
  int i;
  double sma,ecc,incl,node,peri,mnan;
  double e,e1,e2,f,f1,f2,m,de;
  double xx,yy,dxx,dyy,rr;
  double r[3],v[3];
  double ci,si,cn,sn,cw,sw;
  double mjd0=56953.7629284724;

  // Elements
  sma=213180.767;
  ecc=0.96909916;
  incl=28.492086;
  node=299.881688;
  peri=143.175527;
  mnan=0.081443;

  // Bisect Kepler's equation
  m=fabs(mnan)*D2R;
  e1=0.0;
  e2=M_PI;
  f1=e1-ecc*sin(e1)-m;
  f2=e2-ecc*sin(e2)-m;
  for (i=0;i<64;i++) {
    e=0.5*(e1+e2);
    f=e-ecc*sin(e)-m;
    if (f*f1>0.0) {
      f1=f;
      e1=e;
    } else {
      f2=f;
      e2=e;
    }
  }
  e*=fabs(mnan)/mnan;

  // In plane positions
  rr=sma*(1.0-ecc*cos(e));
  de=sqrt(GM_EARTH/sma)/rr;
  xx=sma*(cos(e)-ecc);
  yy=sma*sqrt(1.0-ecc*ecc)*sin(e);
  dxx=-sma*de*sin(e);
  dyy=sma*sqrt(1.0-ecc*ecc)*de*cos(e);

  // Angles
  ci=cos(incl*D2R);
  si=sin(incl*D2R);
  cn=cos(node*D2R);
  sn=sin(node*D2R);
  cw=cos(peri*D2R);
  sw=sin(peri*D2R);
  
  // Position and velocity
  r[0]=xx*(cw*cn-sw*sn*ci)+yy*(-sw*cn-cw*sn*ci);
  r[1]=xx*(cw*sn+sw*cn*ci)+yy*(-sw*sn+cw*cn*ci);
  r[2]=xx*sw*si+yy*cw*si;
  v[0]=dxx*(cw*cn-sw*sn*ci)+dyy*(-sw*cn-cw*sn*ci);
  v[1]=dxx*(cw*sn+sw*cn*ci)+dyy*(-sw*sn+cw*cn*ci);
  v[2]=dxx*sw*si+dyy*cw*si;

  printf("%14.8lf %f %f %f %f %f %f\n",mjd0,r[0],r[1],r[2],v[0],v[1],v[2]);

  return 0;
}
