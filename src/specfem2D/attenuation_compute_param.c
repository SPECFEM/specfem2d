
/* See Liu, Anderson & Kanamori (Geophysical Journal of the Royal Astronomical Society, vol. 47, p. 41-58, 1976) for details */

/* cleaned by Dimitri Komatitsch, University of Pau, France, July 2007 */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* useful constants */

#define PI2 (2 * M_PI)

/* prototypes */

void constant_Q2_sub(double f1, double f2, int n, double Q, double *tau_s, double *tau_e);
void nrerror(const char *error_text);
double *dvector(int n);
double **dmatrix(int nr, int nc);
void free_dvector(double *v, int n);
void free_dmatrix(double **m, int nr, int nc);

/* SVD: Singular value decomposition, see e.g. http://en.wikipedia.org/wiki/Singular_value_decomposition */
void dsvdcmp(double **a, int m, int n, double *w, double **v);

void initialize(double f1, double f2, int n, double Q, double *x1, double *x2);
void derivatives(double f1, double f2, int n, double Q, double *x1, double *x2, double *gradient, double **hessian);
void invert(double *x, double *b, double **A, int n);

/* This is called from file "attenuation_model.f90" */
void
FC_FUNC_(attenuation_compute_param,ATTENUATION_COMPUTE_PARAM)(int *nmech_in,
             double *Q1_in, double *Q2_in,
             double *f1_in, double *f2_in,
             double *tau_sigma_nu1, double *tau_sigma_nu2,
             double *tau_epsilon_nu1, double *tau_epsilon_nu2
             )
{
  int             n;
  double          target_Q1, target_Q2;
  double          f1, f2;

  /* We get the arguments passed from Fortran by address */
  target_Q1 = *Q1_in; /* target value of Q1 */
  target_Q2 = *Q2_in; /* target value of Q2 */
  n = *nmech_in;      /* number of mechanisms */
  f1 = *f1_in;        /* shortest frequency (Hz) */
  f2 = *f2_in;        /* highest frequency (Hz) */

  if (f2 < f1) {
    printf("T2 > T1\n");
    exit(1);
  }

  if (target_Q1 <= -0.0001) {
    printf("Q1 cannot be negative\n");
    exit(1);
  }

  if (target_Q2 <= -0.0001) {
    printf("Q2 cannot be negative\n");
    exit(1);
  }

  if (n < 1) {
    printf("n < 1\n");
    exit(1);
  }

  /* no need to compute these parameters if there is no attenuation; it could lead to a division by zero in the code */
  if (target_Q1 > 0.00001) {
    /* Q1 dilatation mode defined in Carcione's papers */
    constant_Q2_sub(f1, f2, n, target_Q1, tau_sigma_nu1, tau_epsilon_nu1);
  }

  /* no need to compute these parameters if there is no attenuation; it could lead to a division by zero in the code */
  if (target_Q2 > 0.00001) {
    /* Q2 shear mode defined in Carcione's papers */
    constant_Q2_sub(f1, f2, n, target_Q2, tau_sigma_nu2, tau_epsilon_nu2);
  }
}

void constant_Q2_sub(double f1, double f2, int n, double Q, double *tau_s, double *tau_e)
{
  int             i,j;
  double         *x1, *x2;
  double         *gradient, **hessian;

  if (f2 < f1) {
    printf("T2 > T1\n");
    exit(1);
  }
  if (Q < 0.0) {
    printf("Q < 0\n");
    exit(1);
  }
  if (n < 1) {
    printf("n < 1\n");
    exit(1);
  }

  x1 = dvector(n);
  x2 = dvector(n);
  gradient = dvector(n);
  hessian = dmatrix(n, n);
  memset(x1, 0, n*sizeof(double));
  memset(x2, 0, n*sizeof(double));
  memset(gradient, 0, n*sizeof(double));
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++)
      memset(hessian[i], 0, n*sizeof(double));
  }

  initialize(f1, f2, n, Q, x1, x2);

  derivatives(f1, f2, n, Q, x1, x2, gradient, hessian);

  invert(x1, gradient, hessian, n);

  free_dvector(gradient, n);
  free_dmatrix(hessian, n, n);

  for (i = 0; i < n; i++) {
    tau_e[i] = x1[i] + x2[i];
  }
  for (i = 0; i < n; i++) {
    tau_s[i] = x2[i];
  }

  free_dvector(x1, n);
  free_dvector(x2, n);
}

void initialize(double f1, double f2, int n, double Q, double *x1, double *x2)
{
int             i;
double          q, omega, *tau_e, *tau_s;
double          exp1, exp2, dexp, expo;

tau_e = dvector(n);
tau_s = dvector(n);
if (n > 1) {
  exp1 = log10(f1);
  exp2 = log10(f2);
  dexp = (exp2 - exp1) / ((double) (n - 1));
  q = 1.0 / ((n - 1.0) * Q);
  for (i = 0, expo = exp1; i < n; i++, expo += dexp) {
    omega = PI2 * pow(10.0, expo);
    tau_s[i] = 1.0 / omega;
    tau_e[i] = tau_s[i] * (1.0 + q) / (1.0 - q);
  }
} else {
  q = 1.0 / Q;
  exp1 = log10(f1);
  exp2 = log10(f2);
  expo=(exp1+exp2)/2.0;
  omega = PI2 * pow(10.0, expo);
  tau_s[0] = 1.0 / omega;
  tau_e[0] = tau_s[0] * (1.0 + q) / (1.0 - q);
}
/*
 * x1 denotes the parameter tau_e - tau_s and x2 denotes the parameter tau_s
 */
for (i = 0; i < n; i++) {
  x1[i] = tau_e[i] - tau_s[i];
  x2[i] = tau_s[i];
}

free_dvector(tau_e, n);
free_dvector(tau_s, n);
}

void derivatives(double f1, double f2, int n, double Q, double *x1, double *x2, double *gradient, double **hessian)
{
int             i, j;
double          exp1, exp2, dexp, expo;
double          f, df, omega;
double         *dadp, *dbdp, *dqdp, d2qdp2;
double          tau_e, tau_s, a, b, Q_omega;

dadp = dvector(n);
dbdp = dvector(n);
dqdp = dvector(n);
exp1 = log10(f1);
exp2 = log10(f2);
dexp = (exp2 - exp1) / 100.0;
for (i = 0; i < n; i++) {
  gradient[i] = 0.0;
  for (j = 0; j < i; j++) {
    hessian[j][i] = 0.0;
    hessian[j][i] = hessian[i][j];
  }
}
for (expo = exp1; expo <= exp2; expo += dexp) {
  f = pow(10.0, expo);
  df = pow(10.0, expo + dexp) - f;
  omega = PI2 * f;
  a = (double) (1 - n);
  b = 0.0;
  for (i = 0; i < n; i++) {
    tau_e = x1[i] + x2[i];
    tau_s = x2[i];
    a += (1.0 + omega * omega * tau_e * tau_s) /
       (1.0 + omega * omega * tau_s * tau_s);
    b += omega * (tau_e - tau_s) /
    (1.0 + omega * omega * tau_s * tau_s);
    dadp[i] = omega * omega * tau_s / (1.0 + omega * omega * tau_s * tau_s);
    dbdp[i] = omega / (1.0 + omega * omega * tau_s * tau_s);
  }
  Q_omega = a / b;
  for (i = 0; i < n; i++) {
    dqdp[i] = (dbdp[i] - (b / a) * dadp[i]) / a;
    gradient[i] += 2.0 * (1.0 / Q_omega - 1.0 / Q) * dqdp[i] * df / (f2 - f1);
    for (j = 0; j <= i; j++) {
      d2qdp2 = -(dadp[i] * dbdp[j] + dbdp[i] * dadp[j]
           - 2.0 * (b / a) * dadp[i] * dadp[j]) / (a * a);
      hessian[i][j] += (2.0 * dqdp[i] * dqdp[j] + 2.0 * (1.0 / Q_omega - 1.0 / Q) * d2qdp2)
        * df / (f2 - f1);
      hessian[j][i] = hessian[i][j];
    }
  }
}
free_dvector(dadp, n);
free_dvector(dbdp, n);
free_dvector(dqdp, n);
}

void invert(double *x, double *b, double **A, int n)
{
int             i, j, k;
double         *W, **V, **A_inverse;

W = dvector(n);
V = dmatrix(n, n);
A_inverse = dmatrix(n, n);
dsvdcmp(A, n, n, W, V);
for (i = 0; i < n; i++)
  for (j = 0; j < n; j++)
    V[i][j] = (1.0 / W[i]) * A[j][i];
for (i = 0; i < n; i++) {
  for (j = 0; j < n; j++) {
    A_inverse[i][j] = 0.0;
    for (k = 0; k < n; k++)
      A_inverse[i][j] += A[i][k] * V[k][j];
  }
}
free_dvector(W, n);
free_dmatrix(V, n, n);
for (i = 0; i < n; i++) {
  for (j = 0; j < n; j++) {
    x[i] -= A_inverse[i][j] * b[j];
  }
}
free_dmatrix(A_inverse, n, n);
}

void nrerror(const char *error_text)
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

double *dvector(int n)
{
  double *v;

  v=(double *)malloc(n*sizeof(double));
  if (!v) nrerror("allocation failure in dvector()");
  return v;
}

double **dmatrix(int nr, int nc)
{
  int i;
  double **m;

  m=(double **) malloc(nr*sizeof(double*));
  if (!m) nrerror("allocation failure 1 in dmatrix()");

  for(i=0;i<nr;i++) {
    m[i]=(double *) malloc(nc*sizeof(double));
    if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
  }
  return m;
}

void free_dvector(double *v, int n)
{
  free(v);
}

void free_dmatrix(double **m, int nr, int nc)
{
  int i;

  for(i=0;i<nr;i++) free(m[i]);
  free(m);
}

void dsvdcmp(double **a, int m, int n, double *w, double **v)
{
  int flag,i,its,j,jj,k,l,nm;
  double c,f,h,s,x,y,z;
  double anorm=0.0,g=0.0,scale=0.0;
  double *rv1;

  if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
  rv1=dvector(n);
  for (i=0;i<n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a[k][i]);
      if (scale) {
        for (k=i;k<m;k++) {
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f=a[i][i];
        g = -copysign(sqrt(s),f);
        h=f*g-s;
        a[i][i]=f-g;
        if (i != n-1) {
          for (j=l;j<n;j++) {
            for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
            f=s/h;
            for (k=i;k<m;k++) a[k][j] += f*a[k][i];
          }
        }
        for (k=i;k<m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale*g;
    g=s=scale=0.0;
    if (i < m && i != (n-1)) {
      for (k=l;k<n;k++) scale += fabs(a[i][k]);
      if (scale) {
        for (k=l;k<n;k++) {
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g = -copysign(sqrt(s),f);
        h=f*g-s;
        a[i][l]=f-g;
        for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
        if (i != (m-1)) {
          for (j=l;j<m;j++) {
            for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
            for (k=l;k<n;k++) a[j][k] += s*rv1[k];
          }
        }
        for (k=l;k<n;k++) a[i][k] *= scale;
      }
    }
    anorm=fmax(anorm, fabs(w[i])+fabs(rv1[i]));
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g) {
        for (j=l;j<n;j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<n;j++) {
          for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
          for (k=l;k<n;k++) v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=n-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    if (i < n-1)
      for (j=l;j<n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      if (i != n-1) {
        for (j=l;j<n;j++) {
          for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
          f=(s/a[i][i])*g;
          for (k=i;k<m;k++) a[k][j] += f*a[k][i];
        }
      }
      for (j=i;j<m;j++) a[j][i] *= g;
    } else {
      for (j=i;j<m;j++) a[j][i]=0.0;
    }
    ++a[i][i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its<30;its++) {
      flag=1;
      for (l=k;l>=0;l--) {
        nm=l-1;
        if (fabs(rv1[l])+anorm == anorm) {
          flag=0;
          break;
        }
        if (fabs(w[nm])+anorm == anorm) break;
      }
      if (flag) {
        c=0.0;
        s=1.0;
        for (i=l;i<k;i++) {
          f=s*rv1[i];
          if (fabs(f)+anorm != anorm) {
            g=w[i];
            h=hypot(f,g);
            w[i]=h;
            h=1.0/h;
            c=g*h;
            s=(-f*h);
            for (j=0;j<m;j++) {
              y=a[j][nm];
              z=a[j][i];
              a[j][nm]=y*c+z*s;
              a[j][i]=z*c-y*s;
            }
          }
        }
      }
      z=w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j=0;j<n;j++) v[j][k]=(-v[j][k]);
        }
        break;
      }
      if (its == 60) nrerror("No convergence in 60 SVDCMP iterations");
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=hypot(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+copysign(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=hypot(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y=y*c;
        for (jj=0;jj<n;jj++) {
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=hypot(f,h);
        w[j]=z;
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=(c*g)+(s*y);
        x=(c*y)-(s*g);
        for (jj=0;jj<m;jj++) {
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free_dvector(rv1,n);
}

