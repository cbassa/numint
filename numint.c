#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "force.h"
#include "vecmat.h"
#include "refsys.h"
#include "constants.h"
#include "rk.h"
#include "jpleph.h"

double mjd0;
void *ephem;

void accel(double t,double r[],double drdt[])
{
  int i;
  double mjd;
  double e[3][3],f[3][3],rin[3],vin[3],a[3],r_moon[3],r_sun[3],a_moon[3],a_sun[3],a_srp[3],a_drag[3];
  int err_code;
  double area=1e-6,mass=1000.0,cr=1.5,cd=2.0,illum;

  // Get time
  mjd=mjd0+t/86400.0;

  // Compute transformation
  icrs_to_itrs(mjd,e);
  //  eci_to_ecef(mjd,e);
  icrs_to_eme(mjd,f);

  // Input vectors
  for (i=0;i<3;i++) {
    rin[i]=r[i];
    vin[i]=r[i+3];
  }

  // Compute geopotential
  jgm3_geopotential(rin,e,20,20,a);

  // Perturbations by Sun
  sun_position(mjd,ephem,r_sun);
  third_body(rin,r_sun,GM_SUN,a_sun);

  // Perturbations by Moon
  moon_position(mjd,ephem,r_moon);
  third_body(rin,r_moon,GM_MOON,a_moon);

  // Solar radiation pressure
  //  illum=illumination(rin,r_sun);
  //  srp(rin,r_sun,area,mass,cr,a_srp);

  // Drag
  //  drag(rin,vin,r_sun,f,area,mass,cd,a_drag);

  // Derivatives
  for (i=0;i<3;i++) {
    drdt[i]=r[i+3];
    drdt[i+3]=a[i]+a_moon[i];//+a_sun[i]+a_drag[i]+a_srp[i]*illum;
  }

  return;
}

int main(int argc,char *argv[])
{
  int i,j,n=6,flag=0;
  double r[6],drdt[6],rout[6],rerr[6];
  double t,dt,rr;
  char nams[1018][6];
  double vals[1018];

  // Initialize JPL ephemerides
  ephem=jpl_init_ephemeris("lnxp1600p2200.405",nams,vals);

  // CH5T1 injection vector (J2000)
  mjd0=56953.76292847;
  r[0]=-1804.728420;
  r[1]=6494.526142;
  r[2]=906.894086;
  r[3]=-9.709001;
  r[4]=-0.078195;
  r[5]=-4.590375;
  t=0.0;
  dt=120.0;

  for (i=0;i<8640;i++) {
    accel(t,r,drdt);

    printf("%14.8lf %f %f %f %f %f %f\n",mjd0+t/86400.0,r[0],r[1],r[2],r[3],r[4],r[5]);
    rr=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);

    if (rr<XKMPER+100.0) {
      fprintf(stderr,"Satellite reentered Earth's atmosphere\n");
      break;
    }
    rkck(r,drdt,n,t,dt,rout,rerr,accel);
    
    for (j=0;j<n;j++)
      r[j]=rout[j];

    t+=dt;
  }
  
  jpl_close_ephemeris(ephem);

  return 0;
}
