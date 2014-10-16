double modulo(double x,double y);
double gmst(double mjd);
void nutation(double mjd,double *dpsi,double *deps,double *eps);
void precess(double mjd0,double mjd,double *zeta,double *z,double *theta);
void icrs_to_itrs(double mjd,double a[3][3]);
