typedef struct {
        double r;
        double i;
} complex;

typedef struct {
        float r;
        float i;
} f_complex;


double dot(double *a, double *b);
void add(double *a, double *b, double *c);
void cross(double *a, double *b, double *c);
void subtract(double *a, double *b, double *c);
void mixprod(double a, double *b, double *c);
void vfill(double *vect, double a, double b, double c);

float f_dot(float *a, float *b);
void f_add(float *a, float *b, float *c);
void f_cross(float *a, float *b, float *c);
void f_subtract(float *a, float *b, float *c);
void f_mixprod(float a, float *b, float *c);
void f_vfill(float *vect, float a, float b, float c);

complex cdot(complex *a, complex *b);
complex cprod(complex a, complex b);
complex csum(complex a, complex b);
complex cdiff(complex a, complex b);
complex cquot(complex a, complex b);

void cadd(complex *a, complex *b, complex *c);
void ccross(complex *a, complex *b, complex *c);
void csubtract(complex *a, complex *b, complex *c);
void cmixprod(complex a, complex *b, complex *c);
void cvfill(complex *vect, double x1, double x2, double y1, 
            double y2, double z1, double z2);
void cfill(complex *scal, double x1, double x2);

f_complex f_cdot(f_complex *a, f_complex *b);
f_complex f_cprod(f_complex a, f_complex b);
f_complex f_csum(f_complex a, f_complex b);
f_complex f_cdiff(f_complex a, f_complex b);
f_complex f_cquot(f_complex a, f_complex b);
void f_cadd(f_complex *a, f_complex *b, f_complex *c);
void f_ccross(f_complex *a, f_complex *b, f_complex *c);
void f_csubtract(f_complex *a, f_complex *b, f_complex *c);
void f_cmixprod(f_complex a, f_complex *b, f_complex *c);
void f_cvfill(f_complex *vect, float x1, float x2, float y1, 
              float y2, float z1, float z2);
void f_cvfillc(f_complex *vect, f_complex x, f_complex y, f_complex z);
void f_cfill(f_complex *scal, float x1, float x2);
