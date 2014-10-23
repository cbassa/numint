double modulo(double x,double y);
double gmst(double mjd);
void nutation(double mjd,double *dpsi,double *deps,double *eps);
void precess(double mjd0,double mjd,double *zeta,double *z,double *theta);
void icrs_to_itrs(double mjd,double a[3][3]);
void eci_to_ecef(double mjd,double a[3][3]);
void sun_position(double mjd,double r_sun[3]);
void moon_position(double mjd,double r_sun[3]);
