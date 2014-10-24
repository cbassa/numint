#include <stdio.h>
#include <math.h>

// Set identity matrix
void identity_matrix(double a[3][3])
{
  int i,j;

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      if (i==j)
	a[i][j]=1.0;
      else
	a[i][j]=0.0;
    }
  }

  return;
}

// Rotate around x-axis
void rotate_x(double phi, double a[3][3])
{
  double s,c,a10,a11,a12,a20,a21,a22;

  s=sin(phi);
  c=cos(phi);

  a10=c*a[1][0]+s*a[2][0];
  a11=c*a[1][1]+s*a[2][1];
  a12=c*a[1][2]+s*a[2][2];
  a20=-s*a[1][0]+c*a[2][0];
  a21=-s*a[1][1]+c*a[2][1];
  a22=-s*a[1][2]+c*a[2][2];
  
  a[1][0]=a10;
  a[1][1]=a11;
  a[1][2]=a12;
  a[2][0]=a20;
  a[2][1]=a21;
  a[2][2]=a22;
  
  return;
}

// Rotate around y-axis
void rotate_y(double phi, double a[3][3])
{
  double s,c,a00,a01,a02,a20,a21,a22;

  s=sin(phi);
  c=cos(phi);

  a00=c*a[0][0]-s*a[2][0];
  a01=c*a[0][1]-s*a[2][1];
  a02=c*a[0][2]-s*a[2][2];
  a20=s*a[0][0]+c*a[2][0];
  a21=s*a[0][1]+c*a[2][1];
  a22=s*a[0][2]+c*a[2][2];
  
  a[0][0]=a00;
  a[0][1]=a01;
  a[0][2]=a02;
  a[2][0]=a20;
  a[2][1]=a21;
  a[2][2]=a22;
  
  return;
}

// Rotate around z-axis
void rotate_z(double phi, double a[3][3])
{
  double s,c,a00,a01,a02,a10,a11,a12;

  s=sin(phi);
  c=cos(phi);

  a00=c*a[0][0]+s*a[1][0];
  a01=c*a[0][1]+s*a[1][1];
  a02=c*a[0][2]+s*a[1][2];
  a10=-s*a[0][0]+c*a[1][0];
  a11=-s*a[0][1]+c*a[1][1];
  a12=-s*a[0][2]+c*a[1][2];

  a[0][0]=a00;
  a[0][1]=a01;
  a[0][2]=a02;
  a[1][0]=a10;
  a[1][1]=a11;
  a[1][2]=a12;

  return;
}

// Matrix multiply
void matrix_multiply(double a[3][3],double b[3][3],double c[3][3])
{
  int i,j,k;

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      c[i][j]=0.0;
      for (k=0;k<3;k++) 
	c[i][j]+=a[i][k]*b[k][j];
    }
  }

  return;
}

// Vector multiply
void vector_multiply_in_place(double a[3][3],double b[3])
{
  int i,j,k;
  double c[3];

  for (i=0;i<3;i++) {
    c[i]=0.0;
    for (j=0;j<3;j++) 
      c[i]+=a[i][j]*b[j];
  }
  for (i=0;i<3;i++)
    b[i]=c[i];

  return;
}

// Vector multiply
void vector_multiply(double a[3][3],double b[3],double c[3])
{
  int i,j,k;

  for (i=0;i<3;i++) {
    c[i]=0.0;
    for (j=0;j<3;j++) 
      c[i]+=a[i][j]*b[j];
  }

  return;
}

// Transpose
void matrix_transpose(double a[3][3],double b[3][3])
{
  int i,j;

  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      b[i][j]=a[j][i];

  return;
}

// Dot product
double dot_product(double a[3],double b[3])
{
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}  

// Magnitude
double magnitude(double a[3])
{
  return sqrt(dot_product(a,a));
}

// Cross product
void cross_product(double a[3],double b[3],double c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];

  return;
}
