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
  //icrs_to_itrs(mjd,e);
  eci_to_ecef(mjd,e);
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
  illum=illumination(rin,r_sun);
  srp(rin,r_sun,area,mass,cr,a_srp);

  // Drag
  drag(rin,vin,r_sun,f,area,mass,cd,a_drag);

  // Derivatives
  for (i=0;i<3;i++) {
    drdt[i]=r[i+3];
    drdt[i+3]=a[i]+a_moon[i]+a_sun[i];//+a_drag[i]+a_srp[i]*illum;
  }

  return;
}

int main(int argc,char *argv[])
{
  int i,j,n=6;
  double r[6],drdt[6],rout[6],rerr[6];
  double t,dt,rr,rv,ra,r_moon[3],rm;
  char nams[1018][6];
  double vals[1018];

  // Initialize JPL ephemerides
  ephem=jpl_init_ephemeris("lnxp1600p2200.405",nams,vals);

  /*
  // ISS
  mjd0=56946.75000;
  r[0]=5314.1173;
  r[1]=4227.4337;
  r[2]=-222.6590;
  r[3]=-3.11089;
  r[4]=3.60157;
  r[5]=-6.00570;
  t=0.0;
  dt=10.0;
  */
  /*
  // Chang'e 5-T1
  mjd0=56953.78125;
  r[0]=-7510.8025;
  r[1]=+3990.1924;
  r[2]=-2728.2834;
  r[3]=-7.03188;
  r[4]=-4.09134;
  r[5]=-4.57691;
  t=0.0;
  dt=120.0;

  // Apollo 13
  mjd0=40690.36339931;
  r[0]=-226427.771500;
  r[1]=+243502.800515;
  r[2]=+125141.839099;
  r[3]=-0.717752804;
  r[4]=+0.554568851;
  r[5]=+0.257251614;
  t=0.0;
  dt=120.0;


  // Superbird 1A
  mjd0=56293.0;
  r[0]=18585.8274;
  r[1]=34160.0061;
  r[2]=21740.9971;
  r[3]=-1.11886;
  r[4]=-0.45521;
  r[5]=+1.21542;
  t=0.0;
  dt=10.0;

  // COSMOS 2342
  mjd0=56293.0;
  r[0]=18585.9274;
  r[1]=34160.0061;
  r[2]=21740.9971;
  r[3]=-1.11886;
  r[4]=-0.45521;
  r[5]=+1.21542;
  t=0.0;
  dt=120.0;
  */

  // Chang'E 5-T1 spacetrack elset
  mjd0=56954.00000;
  r[0]=-40506.6555;
  r[1]=-60908.2891;
  r[2]=-35723.7912;
  r[3]=-0.67359;
  r[4]=-2.53569;
  r[5]=-0.99681;
  t=0.0;
  dt=120.0;

  for (i=0;i<180000;i++) {
    accel(t,r,drdt);

    moon_position(mjd0+t/86400.0,ephem,r_moon);
    for (j=0;j<3;j++)
      r_moon[j]-=r[j];
    rm=magnitude(r_moon);

    rr=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    rv=sqrt(r[3]*r[3]+r[4]*r[4]+r[5]*r[5]);
    ra=sqrt(drdt[3]*drdt[3]+drdt[4]*drdt[4]+drdt[5]*drdt[5]);
    printf("%14.8lf %f %f %f %f %f %f %f %f %f\n",mjd0+t/86400.0,r[0],r[1],r[2],r[3],r[4],r[5],rr,rv,rm);

    if (rr<XKMPER+100.0) {
      fprintf(stderr,"Satellite crashed\n");
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
