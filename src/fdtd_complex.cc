#include <stdio.h>
#include "fdtd_complex.hh"

double dot(double *a, double *b)
{
  double c;
  int i;
  c=0.;
  for(i=0;i<3;i++)
    c+= a[i]*b[i];
  return(c);
}

void add(double *a, double *b, double *c)
{
  int i;
  for(i=0;i<3;i++)
    c[i]=a[i]+b[i];
  return;
}

void cross(double *a, double *b, double *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  return;
}

void subtract(double *a, double *b, double *c)
{
  int i;
  for(i=0;i<3;i++)
    c[i]=a[i]-b[i];
  return;
}

void mixprod(double a, double *b, double *c)
{
  int i;
  for(i=0;i<3;i++)
    c[i] = a*b[i];
  return;
}

void vfill(double *vect, double x, double y, double z)
{
  vect[0]=x;
  vect[1]=y;
  vect[2]=z;
}

float f_dot(float *a, float *b)
{
  float c;
  int i;
  c=0.;
  for(i=0;i<3;i++)
    c+= a[i]*b[i];
  return(c);
}

void f_add(float *a, float *b, float *c)
{
  int i;
  for(i=0;i<3;i++)
    c[i]=a[i]+b[i];
  return;
}

void f_cross(float *a, float *b, float *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  return;
}

void f_subtract(float *a, float *b, float *c)
{
  int i;
  for(i=0;i<3;i++)
    c[i]=a[i]-b[i];
  return;
}

void f_mixprod(float a, float *b, float *c)
{
  int i;
  for(i=0;i<3;i++)
    c[i] = a*b[i];
  return;
}

void f_vfill(float *vect, float x, float y, float z)
{
  vect[0]=x;
  vect[1]=y;
  vect[2]=z;
}

complex cdot(complex *a, complex *b)
{
  complex c;
  int i;
  c=cprod(*a,*b);
  for(i=1;i<3;i++)
    c=csum(c,cprod(*(a+i),*(b+i)));
  return(c);
}

complex cprod(complex a, complex b)
{
  complex c;
  c.r=a.r*b.r-a.i*b.i;
  c.i=a.r*b.i+a.i*b.r;
  return(c);
}

complex csum(complex a, complex b)
{
  complex c;
  c.r=a.r+b.r;
  c.i=a.i+b.i;
  return(c);
}

complex cdiff(complex a, complex b)
{
  complex c;
  c.r=a.r-b.r;
  c.i=a.i-b.i;
  return(c);
}

void cvfill(complex *vect, double x1, double x2, double y1, 
            double y2, double z1, double z2)
{
  vect->r=x1;
  vect->i=x2;
  (vect+1)->r=y1;
  (vect+1)->i=y2;
  (vect+2)->r=z1;
  (vect+2)->i=z2;
  return;
}

void cfill(complex *scal, double x1, double x2)
{
  scal->r=x1;
  scal->i=x2;
  return;
}

void ccross(complex *a, complex *b, complex *c)
{
  (*c)= cdiff(cprod(*(a+1),*(b+2)),cprod(*(a+2),*(b+1)));
  *(c+1)= cdiff(cprod(*(a+2),*(b)),cprod(*(a),*(b+2)));
  *(c+2)= cdiff(cprod(*(a),*(b+1)),cprod(*(a+1),*(b)));
  return;
}

void cadd(complex *a, complex *b, complex *c)
{
  int i;
  for(i=0;i<3;++i)
    *(c+i)= csum(*(a+i),*(b+i));
  return;
}

void csubtract(complex *a, complex *b, complex *c)
{
  int i;
  for(i=0;i<3;++i)
    *(c+i)= cdiff(*(a+i),*(b+i));
  return;
}

complex cquot(complex a, complex b)
{
  complex c;
  double mod;
  mod= b.r*b.r+b.i*b.i;
  if(mod==0.)printf("dividing by zero ! \n");
  c.r= (a.r*b.r+a.i*b.i)/mod;
  c.i= (a.i*b.r-a.r*b.i)/mod;
  return(c);
}

void cmixprod(complex a, complex *b, complex *c)
{
  int i;
  for(i=0;i<3;++i)
    *(c+i)= cprod(a,*(b+i));
  return;
}

f_complex f_cdot(f_complex *a, f_complex *b)
{
  f_complex c;
  int i;
  c=f_cprod(*a,*b);
  for(i=1;i<3;i++)
    c=f_csum(c,f_cprod(*(a+i),*(b+i)));
  return(c);
}
f_complex f_cprod(f_complex a, f_complex b)
{
  f_complex c;
  c.r=a.r*b.r-a.i*b.i;
  c.i=a.r*b.i+a.i*b.r;
  return(c);
}

f_complex f_csum(f_complex a, f_complex b)
{
  f_complex c;
  c.r=a.r+b.r;
  c.i=a.i+b.i;
  return(c);
}

f_complex f_cdiff(f_complex a, f_complex b)
{
  f_complex c;
  c.r=a.r-b.r;
  c.i=a.i-b.i;
  return(c);
}

void f_cvfill(f_complex *vect, float x1, float x2, float y1, 
              float y2, float z1, float z2)
{
  vect->r=x1;
  vect->i=x2;
  (vect+1)->r=y1;
  (vect+1)->i=y2;
  (vect+2)->r=z1;
  (vect+2)->i=z2;
  return;
}

void f_cvfillc(f_complex *vect, f_complex x, f_complex y, f_complex z)
{
  vect[0]=x;
  vect[1]=y;
  vect[2]=z;
  return;
}

void f_cfill(f_complex *scal, float x1, float x2)
{
  scal->r=x1;
  scal->i=x2;
  return;
}

void f_ccross(f_complex *a, f_complex *b, f_complex *c)
{
  (*c)= f_cdiff(f_cprod(*(a+1),*(b+2)),f_cprod(*(a+2),*(b+1)));
  *(c+1)= f_cdiff(f_cprod(*(a+2),*(b)),f_cprod(*(a),*(b+2)));
  *(c+2)= f_cdiff(f_cprod(*(a),*(b+1)),f_cprod(*(a+1),*(b)));
  return;
}

void f_cadd(f_complex *a, f_complex *b, f_complex *c)
{
  int i;
  for(i=0;i<3;++i)
    *(c+i)= f_csum(*(a+i),*(b+i));
  return;
}

void f_csubtract(f_complex *a, f_complex *b, f_complex *c)
{
  int i;
  for(i=0;i<3;++i)
    *(c+i)= f_cdiff(*(a+i),*(b+i));
  return;
}

f_complex f_cquot(f_complex a, f_complex b)
{
  f_complex c;
  float mod;
  mod= b.r*b.r+b.i*b.i;
  if(mod==0.)printf("dividing by zero ! \n");
  c.r= (a.r*b.r+a.i*b.i)/mod;
  c.i= (a.i*b.r-a.r*b.i)/mod;
  return(c);
}

void f_cmixprod(f_complex a, f_complex *b, f_complex *c)
{
  int i;
  for(i=0;i<3;++i)
    *(c+i)= f_cprod(a,*(b+i));
  return;
}
	
