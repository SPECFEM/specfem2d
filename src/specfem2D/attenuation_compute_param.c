
/* See Liu, Anderson & Kanamori (Geophysical Journal of the Royal Astronomical Society, vol. 47, p. 41-58, 1976) for details */

/* cleaned by Dimitri Komatitsch, University of Pau, France, July 2007 */

#include "config.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include <sgtty.h>
#include <signal.h>
#include <stdlib.h>

/* useful constants */

#define PI2 (2 * M_PI)

/* It is called in "attenuation_model.f90". */
void
FC_FUNC_(attenuation_compute_param,ATTENUATION_COMPUTE_PARAM)(int *nmech_in,
             double *Q1_in, double *Q2_in,
             double *f1_in, double *f2_in,
             double *tau_sigma_nu1, double *tau_sigma_nu2,
             double *tau_epsilon_nu1, double *tau_epsilon_nu2
             )
{
  int             n, i, j, plot, nu;
  double          Q_value, target_Q1, target_Q2;
  double          f1, f2, Q, om0, Omega;
  double          a, b;
  double         *tau_s, *tau_e;
  double         *dvector();
  void            constant_Q2_sub();
  void            free_dvector();


  /* We get the arguments passed in fortran by adress. */
  target_Q1 = *Q1_in; /* target value of Q1 */
  target_Q2 = *Q2_in; /* target value of Q2 */
  n = *nmech_in;      /* number of mechanisms */
  f1 = *f1_in;        /* shortest frequency (Hz) */
  f2 = *f2_in;        /* highest frequency (Hz) */

  if (f2 < f1) {
    printf("T2 > T1\n");
    exit; }

  if (target_Q1 <= -0.0001) {
    printf("Q1 cannot be negative\n");
    exit; }

  if (target_Q2 <= -0.0001) {
    printf("Q2 cannot be negative\n");
    exit; }

  if (n < 1) {
    printf("n < 1\n");
    exit; }

  om0 = PI2 * pow(10.0, 0.5 * (log10(f1) + log10(f2)));

  plot = 0;

/* loop on the Q1 dilatation mode (nu = 1) and Q2 shear mode (nu = 2) defined in Carcione's papers */
  for (nu = 1; nu <= 2; nu++) {

/* assign Q1 or Q2 to generic variable Q_value which is used for the calculations */
    if (nu == 1) { Q_value = target_Q1 ; }
    if (nu == 2) { Q_value = target_Q2 ; }

/* no need to compute these parameters if there is no attenuation; it could lead to a division by zero in the code */
    if (Q_value > 0.00001) {

    tau_s = dvector(1, n);
    tau_e = dvector(1, n);

    constant_Q2_sub(f1, f2, n, Q_value, tau_s, tau_e);

/* output in Fortran90 format */
    for (i = 1; i <= n; i++) {
      /* We put the results in tau_sigma_nu to get them in fortran. */
      if ( nu == 1 ) {
        tau_sigma_nu1[i-1] = tau_s[i];
      }
      if ( nu == 2 ) {
        tau_sigma_nu2[i-1] = tau_s[i];
      }

    }

    for (i = 1; i <= n; i++) {
       /* We put the results in tau_epsilon_nu to get them in fortran. */
      if ( nu == 1 ) {
        tau_epsilon_nu1[i-1] = tau_e[i];
      }
      if ( nu == 2 ) {
        tau_epsilon_nu2[i-1] = tau_e[i];
      }

    }

    free_dvector(tau_s, 1, n);
    free_dvector(tau_e, 1, n);
  }
  }
}

#include <stdio.h>

void nrerror(error_text)
char error_text[];
{
  void exit();

  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

double *dvector(nl,nh)
int nl,nh;
{
  double *v;

  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl;
}

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
  int i;
  double **m;

  m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
  free((char*) (v+nl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
  int i;

  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

#include <math.h>

static double at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static double maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
  (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void dsvdcmp(a,m,n,w,v)
double **a,*w,**v;
int m,n;
{
  int flag,i,its,j,jj,k,l,nm;
  double c,f,h,s,x,y,z;
  double anorm=0.0,g=0.0,scale=0.0;
  double *rv1,*dvector();
  void nrerror(),free_dvector();

  if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
  rv1=dvector(1,n);
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
        for (k=i;k<=m;k++) {
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f=a[i][i];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][i]=f-g;
        if (i != n) {
          for (j=l;j<=n;j++) {
            for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
            f=s/h;
            for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
          }
        }
        for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
        for (k=l;k<=n;k++) {
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][l]=f-g;
        for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
        if (i != m) {
          for (j=l;j<=m;j++) {
            for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
            for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
          }
        }
        for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
        for (j=l;j<=n;j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<=n;j++) {
          for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
          for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=n;i>=1;i--) {
    l=i+1;
    g=w[i];
    if (i < n)
      for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      if (i != n) {
        for (j=l;j<=n;j++) {
          for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
          f=(s/a[i][i])*g;
          for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
        }
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else {
      for (j=i;j<=m;j++) a[j][i]=0.0;
    }
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
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
        for (i=l;i<=k;i++) {
          f=s*rv1[i];
          if (fabs(f)+anorm != anorm) {
            g=w[i];
            h=PYTHAG(f,g);
            w[i]=h;
            h=1.0/h;
            c=g*h;
            s=(-f*h);
            for (j=1;j<=m;j++) {
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
          for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
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
      g=PYTHAG(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=PYTHAG(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y=y*c;
        for (jj=1;jj<=n;jj++) {
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=PYTHAG(f,h);
        w[j]=z;
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=(c*g)+(s*y);
        x=(c*y)-(s*g);
        for (jj=1;jj<=m;jj++) {
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
  free_dvector(rv1,1,n);
}

#undef SIGN
#undef MAX
#undef PYTHAG
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include <sgtty.h>
#include <signal.h>
#include <stdlib.h>

void constant_Q2_sub(f1, f2, n, Q, tau_s, tau_e)

  int             n;
  double          f1, f2, Q;
  double         *tau_s, *tau_e;
{
  int             i,j;
  double         *x1, *x2;
  double         *gradient, **hessian;
  double         *dvector(), **dmatrix();
  void            derivatives();
  void            initialize(), invert();
  void            free_dvector(), free_dmatrix();

  if (f2 < f1) {
    printf("T2 > T1\n");
    exit;
  }
  if (Q < 0.0) {
    printf("Q < 0\n");
    exit;
  }
  if (n < 1) {
    printf("n < 1\n");
    exit;
  }

  x1 = dvector(1, n);
  x2 = dvector(1, n);
  gradient = dvector(1, n);
  hessian = dmatrix(1, n, 1, n);
  for(i=1;i<=n;i++) {
    x1[i]=0.0;
    x2[i]=0.0;
    gradient[i]=0.0;
    for(j=1;j<=n;j++) hessian[i][j]=0.0;
  }

  initialize(f1, f2, n, Q, x1, x2);

  derivatives(f1, f2, n, Q, x1, x2, gradient, hessian);

  invert(x1, gradient, hessian, n);

  free_dvector(gradient, 1, n);
  free_dmatrix(hessian, 1, n, 1, n);

  for (i = 1; i <= n; i++) {
          tau_e[i]=x1[i] + x2[i];
  }
  for (i = 1; i <= n; i++) {
          tau_s[i]=x2[i];
  }

  free_dvector(x1, 1, n);
  free_dvector(x2, 1, n);

}

void            initialize(f1, f2, n, Q, x1, x2)
  int             n;
  double          f1, f2, Q, *x1, *x2;
{
int             i;
double          q, omega, *tau_e, *tau_s;
double          exp1, exp2, dexp, expo;
double         *dvector();
void            free_dvector();

tau_e = dvector(1, n);
tau_s = dvector(1, n);
if (n > 1) {
  exp1 = log10(f1);
  exp2 = log10(f2);
  dexp = (exp2 - exp1) / ((double) (n - 1));
  q = 1.0 / ((n - 1.0) * Q);
  for (i = 1, expo = exp1; i <= n; i++, expo += dexp) {
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
  tau_s[1] = 1.0 / omega;
  tau_e[1] = tau_s[1] * (1.0 + q) / (1.0 - q);
}
/*
 * x1 denotes the parameter tau_e - tau_s and x2 denotes the parameter tau_s
 */
for (i = 1; i <= n; i++) {
  x1[i] = tau_e[i] - tau_s[i];
  x2[i] = tau_s[i];
}

free_dvector(tau_e, 1, n);
free_dvector(tau_s, 1, n);
}

void            derivatives(f1, f2, n, Q, x1, x2, gradient, hessian)
  int             n;
  double          f1, f2, Q, *x1, *x2;
  double         *gradient, **hessian;
{
int             i, j;
double          exp1, exp2, dexp, expo;
double          f, df, omega;
double         *dadp, *dbdp, *dqdp, d2qdp2;
double          tau_e, tau_s, a, b, Q_omega;
double         *dvector();
void            free_dvector();

dadp = dvector(1, n);
dbdp = dvector(1, n);
dqdp = dvector(1, n);
exp1 = log10(f1);
exp2 = log10(f2);
dexp = (exp2 - exp1) / 100.0;
for (i = 1; i <= n; i++) {
  gradient[i] = 0.0;
  for (j = 1; j <= i; j++) {
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
  for (i = 1; i <= n; i++) {
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
  for (i = 1; i <= n; i++) {
    dqdp[i] = (dbdp[i] - (b / a) * dadp[i]) / a;
    gradient[i] += 2.0 * (1.0 / Q_omega - 1.0 / Q) * dqdp[i] * df / (f2 - f1);
    for (j = 1; j <= i; j++) {
      d2qdp2 = -(dadp[i] * dbdp[j] + dbdp[i] * dadp[j]
           - 2.0 * (b / a) * dadp[i] * dadp[j]) / (a * a);
      hessian[i][j] += (2.0 * dqdp[i] * dqdp[j] + 2.0 * (1.0 / Q_omega - 1.0 / Q) * d2qdp2)
        * df / (f2 - f1);
      hessian[j][i] = hessian[i][j];
    }
  }
}
free_dvector(dadp, 1, n);
free_dvector(dbdp, 1, n);
free_dvector(dqdp, 1, n);
}

void            invert(x, b, A, n)
  int             n;
  double         *x;
  double         *b, **A;
{
int             i, j, k;
double         *dvector(), **dmatrix();
double         *xp, *W, **V, **A_inverse;
void            free_dvector(), free_dmatrix(), dsvdcmp();

xp = dvector(1, n);
W = dvector(1, n);
V = dmatrix(1, n, 1, n);
A_inverse = dmatrix(1, n, 1, n);
dsvdcmp(A, n, n, W, V);
for (i = 1; i <= n; i++)
  for (j = 1; j <= n; j++)
    V[i][j] = (1.0 / W[i]) * A[j][i];
for (i = 1; i <= n; i++) {
  for (j = 1; j <= n; j++) {
    A_inverse[i][j] = 0.0;
    for (k = 1; k <= n; k++)
      A_inverse[i][j] += A[i][k] * V[k][j];
  }
}
free_dvector(W, 1, n);
free_dmatrix(V, 1, n, 1, n);
for (i = 1; i <= n; i++) {
  xp[i] = x[i];
  for (j = 1; j <= n; j++) {
    xp[i] -= A_inverse[i][j] * b[j];
  }
  x[i] = xp[i];
}
free_dvector(xp, 1, n);
free_dmatrix(A_inverse, 1, n, 1, n);
}
