#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "sgdp4h.h"
#include "satutl.h"
#include <getopt.h>
#include "refsys.h"
#include "jpleph.h"
#include "vecmat.h"

#define LIM 128
#define XKE 0.07436680 // Guassian Gravitational Constant
#define XKMPER 6378.135
#define AE 1.0
#define XMNPDA 1440.0
#define R2D 180.0/M_PI
#define D2R M_PI/180.0

extern double SGDP4_jd0;

// Present nfd
void nfd_now(char *s)
{
  time_t rawtime;
  struct tm *ptm;

  // Get UTC time
  time(&rawtime);
  ptm=gmtime(&rawtime);
    
  sprintf(s,"%04d-%02d-%02dT%02d:%02d:%02d",ptm->tm_year+1900,ptm->tm_mon+1,ptm->tm_mday,ptm->tm_hour,ptm->tm_min,ptm->tm_sec);
  
  return;
}

// Compute Julian Day from Date
double date2mjd(int year,int month,double day)
{
  int a,b;
  double jd;

  if (month<3) {
    year--;
    month+=12;
  }

  a=floor(year/100.);
  b=2.-a+floor(a/4.);

  if (year<1582) b=0;
  if (year==1582 && month<10) b=0;
  if (year==1852 && month==10 && day<=4) b=0;

  jd=floor(365.25*(year+4716))+floor(30.6001*(month+1))+day+b-1524.5;

  return jd-2400000.5;
}

// nfd2mjd
double nfd2mjd(char *date)
{
  int year,month,day,hour,min,sec;
  double mjd,dday;

  sscanf(date,"%04d-%02d-%02dT%02d:%02d:%02d",&year,&month,&day,&hour,&min,&sec);
  dday=day+hour/24.0+min/1440.0+sec/86400.0;

  mjd=date2mjd(year,month,dday);

  return mjd;
}

void usage(void)
{

  return;
}

int main(int argc,char *argv[])
{
  int imode,arg,satno=0,useepoch=0,format=0;
  FILE *file;
  char *env;
  char tlefile[LIM],nfd[32];
  double mjd;
  xyz_t r,v;
  orbit_t orb;
  double rr[3],vv[3],e[3][3],et[3][3];

  // Get environment variable
  env=getenv("ST_TLEDIR");
  sprintf(tlefile,"%s/classfd.tle",env);

  // Set date
  nfd_now(nfd);
  mjd=nfd2mjd(nfd);

  // Decode options
  while ((arg=getopt(argc,argv,"c:i:t:m:hef"))!=-1) {
    switch (arg) {

    case 't':
      strcpy(nfd,optarg);
      mjd=nfd2mjd(nfd);
      break;
    
    case 'm':
      mjd=atof(optarg);
      break;

    case 'e':
      useepoch=1;
      break;

    case 'f':
      format=1;
      break;

    case 'c':
      strcpy(tlefile,optarg);
      break;

    case 'i':
      satno=atoi(optarg);
      break;

    case 'h':
      usage();
      return 0;
      break;

    default:
      usage();
      return 0;
    }
  }

  // Open file
  file=fopen(tlefile,"r");
  while (read_twoline(file,satno,&orb)==0) {
    // Propagate
    imode=init_sgdp4(&orb);

    // Use epoch instead of user supplied date
    if (useepoch==1)
      mjd=SGDP4_jd0-2400000.5;

    // Compute position and velocity
    imode=satpos_xyz(mjd+2400000.5,&r,&v);

    // Output
    //    printf("%14.8lf %f %f %f %f %f %f TEME\n",mjd,r.x,r.y,r.z,v.x,v.y,v.z);

    // To vectors
    rr[0]=r.x;
    rr[1]=r.y;
    rr[2]=r.z;
    vv[0]=v.x;
    vv[1]=v.y;
    vv[2]=v.z;

    // Matrices
    icrs_to_teme(mjd,e);
    matrix_transpose(e,et);

    // Transform
    //    vector_multiply_in_place(et,rr);
    //    vector_multiply_in_place(et,vv);
    
    // Output J2000
    if (format==0)
      printf("%14.8lf %f %f %f %f %f %f J2000\n",mjd,rr[0],rr[1],rr[2],vv[0],vv[1],vv[2]);
    else
      printf("mjd0=%14.8lf;\nr[0]=%f;\nr[1]=%f;\nr[2]=%f;\nr[3]=%f;\nr[4]=%f;\nr[5]=%f;\n",mjd,rr[0],rr[1],rr[2],vv[0],vv[1],vv[2]);
  }
  fclose(file);

  return 0;
}
