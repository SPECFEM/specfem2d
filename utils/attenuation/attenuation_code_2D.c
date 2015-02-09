
/* See Liu, Anderson & Kanamori (GJRAS, 47, 41-58, 1976) for details */

/* cleaned by Dimitri Komatitsch, University of Pau, France, July 2007 */

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include <sgtty.h>
#include <signal.h>
#include <stdlib.h>

/* useful constants */

#define PI 3.14159265358979
#define PI2 6.28318530717958

main (argc,argv)
int argc; char **argv;
{
  int             xmgr, n, i, j, plot, nu;
  double          Q_sigma, target_QKappa, target_Qmu;
  double          f1, f2, Q, om0, Omega;
  double          a, b;
  double          kappa, mu, kappa0, mu0, kappaR, muR;
  double         *tau_sigma, *tau_epsilon;
  double         *dvector();
  void            constant_Q2_sub(),plot_modulus();
  void            free_dvector();

  printf("target value of QKappa: ");
  scanf("%lf",&target_QKappa);
  printf("%lf\n",target_QKappa);

  printf("target value of Qmu: ");
  scanf("%lf",&target_Qmu);
  printf("%lf\n",target_Qmu);

  printf("shortest frequency (Hz): ");
  scanf("%lf",&f1);
  printf("%lf\n",f1);

  printf("highest frequency (Hz): ");
  scanf("%lf",&f2);
  printf("%lf\n",f2);

  printf("number of mechanisms: ");
  scanf("%d",&n);
  printf("%d\n",n);

/*  DK DK  printf("1 = use xmgr  0 = do not use xmgr: "); */
/*  scanf("%d",&xmgr);  */
  xmgr = 0;

  if (f2 < f1) {
    printf("T2 > T1\n");
    exit; }

  if (target_QKappa <= 0.0001) {
    printf("QKappa cannot be negative or null\n");
    exit; }

  if (target_Qmu <= 0.0001) {
    printf("Qmu cannot be negative or null\n");
    exit; }

  if (n < 1) {
    printf("n < 1\n");
    exit; }

  om0 = PI2 * pow(10.0, 0.5 * (log10(f1) + log10(f2)));

  printf("\n! put this in file constants.h\n\n");

  printf("! number of standard linear solids for attenuation\n");
  printf("  integer, parameter :: N_SLS = %d\n\n",n);

  printf("! put this in file attenuation_model.f90\n\n");

  printf("! frequency range: %lf Hz - %lf Hz\n", f1 , f2);
  printf("! central frequency in log scale in Hz = %20.15f\n",om0 / PI2);

  printf("! target constant attenuation factor QKappa = %20.10lf\n", target_QKappa);
  printf("! target constant attenuation factor Qmu = %20.10lf\n\n", target_Qmu);

  printf("! tau_sigmaigma evenly spaced in log frequency, do not depend on value of Q\n\n");

  plot = 0;

/* loop on the QKappa dilatation mode (nu = 1) and Qmu shear mode (nu = 2) */
  for (nu = 1; nu <= 2; nu++) {

/* assign QKappa or Qmu to generic variable Q_sigma which is used for the calculations */
    if (nu == 1) { Q_sigma = target_QKappa ; }
    if (nu == 2) { Q_sigma = target_Qmu ; }

    tau_sigma = dvector(1, n);
    tau_epsilon = dvector(1, n);

    constant_Q2_sub(f1, f2, n, Q_sigma, tau_sigma, tau_epsilon, xmgr);

/* output in Fortran90 format */
    for (i = 1; i <= n; i++) {
      printf("  tau_sigmaigma_nu%d(%1d) = %30.20lfd0\n", nu, i, tau_sigma[i]);
      }
      printf("\n");

    for (i = 1; i <= n; i++) {
      printf("  tau_epsilonpsilon_nu%d(%1d) = %30.20lfd0\n", nu, i, tau_epsilon[i]);
      }
    printf("\n");

    free_dvector(tau_sigma, 1, n);
    free_dvector(tau_epsilon, 1, n);

  }

}

void   plot_modulus(f1, f2, n, m, mR, Q, tau_epsilon, tau_sigma ,xmgr)
        int  n, xmgr;
        double f1, f2, m, mR, Q, *tau_epsilon, *tau_sigma;
{
int             pid, i;
double          exp1, exp2, dexp, expo;
double          f, om, Omega;
double          a, b, m_om, m_prem;
char            strng[180];
int             getpid(), system();
FILE           *fp_v, *fp_q;

pid = getpid();
sprintf(strng, "modulus%1d", pid);
if((fp_v=fopen(strng,"w"))==NULL) {
  puts("cannot open file\n");
  exit;
}
sprintf(strng, "Q%1d", pid);
if((fp_q=fopen(strng,"w"))==NULL) {
  puts("cannot open file\n");
  exit;
}

exp1 = log10(f1) - 2.0;
exp2 = log10(f2) + 2.0;
dexp = (exp2 - exp1) / 100.0;
for (expo = exp1; expo <= exp2; expo += dexp) {
  f = pow(10.0, expo);
  om = PI2 * f;
        a = 1.0;
        b = 0.0;
        for (i = 1; i <= n; i++) {
            a -= om * om * tau_epsilon[i] * (tau_epsilon[i] - tau_sigma[i]) /
                (1.0 + om * om * tau_epsilon[i] * tau_epsilon[i]);
          b += om * (tau_epsilon[i] - tau_sigma[i]) /
             (1.0 + om * om * tau_epsilon[i] * tau_epsilon[i]);
        }
        Omega=a*(sqrt(1.0+b*b/(a*a))-1.0);
        m_om = 2.0*mR* Omega/(b*b);
        m_prem = m * (1.0 + (2.0 / (PI * Q)) * log(om / PI2));
        fprintf(fp_v, "%f %f %f\n", expo, m_om/m, m_prem/m);
  if (om >= PI2 * f1 && om <= PI2 * f2) {
           fprintf(fp_q, "%f %f %f\n", expo, 1.0/atan(b/a), Q);
        }
}
fclose(fp_v);
fclose(fp_q);

/* DK DK call xmgr to plot curves if needed */

if (xmgr == 1) {
  sprintf(strng, "xmgr -nxy Q%1d", pid);
  system(strng);
  sprintf(strng, "xmgr -nxy modulus%1d", pid);
  system(strng);
  sprintf(strng, "rm modulus%1d", pid);
  system(strng);
  sprintf(strng, "rm Q%1d", pid);
  system(strng);
}

}

#include <malloc.h>
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

float *vector(nl,nh)
int nl,nh;
{
  float *v;

  v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl;
}

int *ivector(nl,nh)
int nl,nh;
{
  int *v;

  v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl;
}

double *dvector(nl,nh)
int nl,nh;
{
  double *v;

  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl;
}



float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
  int i;
  float **m;

  m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
    if (!m[i]) nrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
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

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
  int i,**m;

  m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
  if (!m) nrerror("allocation failure 1 in imatrix()");
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
    if (!m[i]) nrerror("allocation failure 2 in imatrix()");
    m[i] -= ncl;
  }
  return m;
}



float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
  int i,j;
  float **m;

  m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*));
  if (!m) nrerror("allocation failure in submatrix()");
  m -= newrl;

  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

  return m;
}



void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
  free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v,nl,nh;
{
  free((char*) (v+nl));
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
  free((char*) (v+nl));
}



void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
{
  int i;

  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
  int i;

  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
  int i;

  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}



void free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
  free((char*) (b+nrl));
}



float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
int nrl,nrh,ncl,nch;
{
  int i,j,nrow,ncol;
  float **m;

  nrow=nrh-nrl+1;
  ncol=nch-ncl+1;
  m = (float **) malloc((unsigned) (nrow)*sizeof(float*));
  if (!m) nrerror("allocation failure in convert_matrix()");
  m -= nrl;
  for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
  return m;
}



void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
  free((char*) (b+nrl));
}

#include <math.h>

#define NMAX 5000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

#define GET_PSUM for (j=1;j<=ndim;j++) { for (i=1,sum=0.0;i<=mpts;i++)\
            sum += p[i][j]; psum[j]=sum;}

void amoeba(p,y,ndim,ftol,funk,nfunk)
float **p,y[],ftol,(*funk)();
int ndim,*nfunk;
{
  int i,j,ilo,ihi,inhi,mpts=ndim+1;
  float ytry,ysave,sum,rtol,amotry(),*psum,*vector();
  void nrerror(),free_vector();

  psum=vector(1,ndim);
  *nfunk=0;
  GET_PSUM
  for (;;) {
    ilo=1;
    ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
    for (i=1;i<=mpts;i++) {
      if (y[i] < y[ilo]) ilo=i;
      if (y[i] > y[ihi]) {
        inhi=ihi;
        ihi=i;
      } else if (y[i] > y[inhi])
        if (i != ihi) inhi=i;
    }
    rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
    if (rtol < ftol) break;
    if (*nfunk >= NMAX) nrerror("Too many iterations in AMOEBA");
    ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,-ALPHA);
    if (ytry <= y[ilo])
      ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,GAMMA);
    else if (ytry >= y[inhi]) {
      ysave=y[ihi];
      ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,BETA);
      if (ytry >= ysave) {
        for (i=1;i<=mpts;i++) {
          if (i != ilo) {
            for (j=1;j<=ndim;j++) {
              psum[j]=0.5*(p[i][j]+p[ilo][j]);
              p[i][j]=psum[j];
            }
            y[i]=(*funk)(psum);
          }
        }
        *nfunk += ndim;
        GET_PSUM
      }
    }
  }
  free_vector(psum,1,ndim);
}

float amotry(p,y,psum,ndim,funk,ihi,nfunk,fac)
float **p,*y,*psum,(*funk)(),fac;
int ndim,ihi,*nfunk;
{
  int j;
  float fac1,fac2,ytry,*ptry,*vector();
  void nrerror(),free_vector();

  ptry=vector(1,ndim);
  fac1=(1.0-fac)/ndim;
  fac2=fac1-fac;
  for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
  ytry=(*funk)(ptry);
  ++(*nfunk);
  if (ytry < y[ihi]) {
    y[ihi]=ytry;
    for (j=1;j<=ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  free_vector(ptry,1,ndim);
  return ytry;
}

#undef ALPHA
#undef BETA
#undef GAMMA
#undef NMAX

void spline(x,y,n,yp1,ypn,y2)
float x[],y[],yp1,ypn,y2[];
int n;
{
  int i,k;
  float p,qn,sig,un,*u,*vector();
  void free_vector();

  u=vector(1,n-1);
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free_vector(u,1,n-1);
}

void splint(xa,ya,y2a,n,x,y)
float xa[],ya[],y2a[],x,*y;
int n;
{
  int klo,khi,k;
  float h,b,a;
  void nrerror();

  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

#define FUNC(x) ((*func)(x))

float trapzd(func,a,b,n)
float a,b;
float (*func)();  /* ANSI: float (*func)(float); */
int n;
{
  float x,tnm,sum,del;
  static float s;
  static int it;
  int j;

  if (n == 1) {
    it=1;
    return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
  } else {
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
    it *= 2;
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}

#include <math.h>

#define EPS 0.5e-5
#define JMAX 20
#define JMAXP JMAX+1
#define K 5

float qromb(func,a,b)
float a,b;
float (*func)();
{
  float ss,dss,trapzd();
  float s[JMAXP+1],h[JMAXP+1];
  int j;
  void polint(),nrerror();

  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=trapzd(func,a,b,j);
    if (j >= K) {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) < EPS*fabs(ss)) return ss;
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j];
  }
  nrerror("Too many steps in routine QROMB");
/* return an error code if needed */
  return -1.;
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#include <math.h>

void polint(xa,ya,n,x,y,dy)
float xa[],ya[],x,*y,*dy;
int n;
{
  int i,m,ns=1;
  float den,dif,dift,ho,hp,w;
  float *c,*d,*vector();
  void nrerror(),free_vector();

  dif=fabs(x-xa[1]);
  c=vector(1,n);
  d=vector(1,n);
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) nrerror("Error in routine POLINT");
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free_vector(d,1,n);
  free_vector(c,1,n);
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(idum)
int *idum;
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;

  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

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

/* useful constants */

#define PI 3.14159265358979
#define PI2 6.28318530717958

void constant_Q2_sub(f1, f2, n, Q, tau_sigma, tau_epsilon, xmgr)

  int             n, xmgr;
  double          f1, f2, Q;
  double         *tau_sigma, *tau_epsilon;
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
          tau_epsilon[i]=x1[i] + x2[i];
  }
  for (i = 1; i <= n; i++) {
          tau_sigma[i]=x2[i];
  }

  free_dvector(x1, 1, n);
  free_dvector(x2, 1, n);

}

void            initialize(f1, f2, n, Q, x1, x2)
  int             n;
  double          f1, f2, Q, *x1, *x2;
{
int             i;
double          q, omega, *tau_epsilon, *tau_sigma;
double          exp1, exp2, dexp, expo;
double         *dvector();
void            free_dvector();

tau_epsilon = dvector(1, n);
tau_sigma = dvector(1, n);
if (n > 1) {
  exp1 = log10(f1);
  exp2 = log10(f2);
  dexp = (exp2 - exp1) / ((double) (n - 1));
  q = 1.0 / ((n - 1.0) * Q);
  for (i = 1, expo = exp1; i <= n; i++, expo += dexp) {
    omega = PI2 * pow(10.0, expo);
    tau_sigma[i] = 1.0 / omega;
    tau_epsilon[i] = tau_sigma[i] * (1.0 + q) / (1.0 - q);
  }
} else {
  q = 1.0 / Q;
  exp1 = log10(f1);
  exp2 = log10(f2);
    expo=(exp1+exp2)/2.0;
  omega = PI2 * pow(10.0, expo);
  tau_sigma[1] = 1.0 / omega;
  tau_epsilon[1] = tau_sigma[1] * (1.0 + q) / (1.0 - q);
}
/*
 * x1 denotes the parameter tau_epsilon - tau_sigma and x2 denotes the parameter tau_sigma
 */
for (i = 1; i <= n; i++) {
  x1[i] = tau_epsilon[i] - tau_sigma[i];
  x2[i] = tau_sigma[i];
}

free_dvector(tau_epsilon, 1, n);
free_dvector(tau_sigma, 1, n);
}

double          penalty(f1, f2, n, Q, x1, x2)
  int             n;
  double          f1, f2, Q, *x1, *x2;
{
int             i;
double          exp1, exp2, dexp, expo;
double          pnlt;
double          f, df, omega;
double          tau_epsilon, tau_sigma, a, b, Q_omega;

exp1 = log10(f1);
exp2 = log10(f2);
dexp = (exp2 - exp1) / 100.0;
pnlt = 0.0;
for (expo = exp1; expo <= exp2; expo += dexp) {
  f = pow(10.0, expo);
  df = pow(10.0, expo + dexp) - f;
  omega = PI2 * f;
  a = (double) (1 - n);
  b = 0.0;
  for (i = 1; i <= n; i++) {
    tau_epsilon = x1[i] + x2[i];
    tau_sigma = x2[i];
    a += (1.0 + omega * omega * tau_epsilon * tau_sigma) /
       (1.0 + omega * omega * tau_sigma * tau_sigma);
    b += omega * (tau_epsilon - tau_sigma) /
       (1.0 + omega * omega * tau_sigma * tau_sigma);
  }
  Q_omega = a / b;
  pnlt += pow(1.0 / Q - 1.0 / Q_omega, 2.0) * df;
}
pnlt /= (f2 - f1);
return pnlt;
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
double          tau_epsilon, tau_sigma, a, b, Q_omega;
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
    tau_epsilon = x1[i] + x2[i];
    tau_sigma = x2[i];
    a += (1.0 + omega * omega * tau_epsilon * tau_sigma) /
       (1.0 + omega * omega * tau_sigma * tau_sigma);
    b += omega * (tau_epsilon - tau_sigma) /
    (1.0 + omega * omega * tau_sigma * tau_sigma);
    dadp[i] = omega * omega * tau_sigma / (1.0 + omega * omega * tau_sigma * tau_sigma);
    dbdp[i] = omega / (1.0 + omega * omega * tau_sigma * tau_sigma);
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
