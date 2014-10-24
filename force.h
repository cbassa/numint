void jgm3_geopotential(double r[3],double e[3][3],int nmax,int mmax,double a[3]);
void srp(double r[3],double s[3],double area,double mass,double cr,double a[3]);
void third_body(double r[3],double s[3],double gm,double a[3]);
double illumination(double r[3],double s[3]);
void drag(double r[3],double v[3],double s[3],double e[3][3],double area,double mass,double cd,double a[3]);
double density(double r[3],double s[3]);
