void identity_matrix(double a[3][3]);
void rotate_x(double phi, double a[3][3]);
void rotate_y(double phi, double a[3][3]);
void rotate_z(double phi, double a[3][3]);
void matrix_multiply(double a[3][3],double b[3][3],double c[3][3]);
void matrix_transpose(double a[3][3],double b[3][3]);
void vector_multiply(double a[3][3],double b[3],double c[3]);
void vector_multiply_in_place(double a[3][3],double b[3]);
double dot_product(double a[3],double b[3]);
double magnitude(double a[3]);