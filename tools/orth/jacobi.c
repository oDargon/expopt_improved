/**************************************************************************/
/*                                                                        */
/*                              J A C O B I                               */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* This subroutine diagonalizes the symmetric matrix f, which is stored   */
/* in lower triangular form. The dimension of f is defined by n. The same */
/* unitary transformation matrix is applied to the matrix u. u is defined */
/* as a set of column vectors of dimension nr, so u have the dimension    */
/* nr by n. For normal usage, u is the unit matrix of dimension n on      */
/* entry, and will in that case contain the eigenvectors to f.            */
/*                                                                        */
/**************************************************************************/
#include <math.h>
void jacobi(double f[], double u[], int n, int nr) {
   int i,j,k;
   int ii,jj,ij;
   int ki,kj;
   double eps,err;
   double fii,fjj,fij;
   double diff,sign;
   double d,t,c,s,tmp;
   eps=1.0e-20;
   err=1.0;
   while (err>eps) {
      err=0.0;
      for (i=1; i<n; i++) {
         for (j=0; j<i; j++) {
            ii=i+i*(i+1)/2;
            jj=j+j*(j+1)/2;
            ij=j+i*(i+1)/2;
            fii=f[ii];
            fjj=f[jj];
            fij=f[ij];
            tmp=fij;
            if(fij<0) tmp=-fij;
            if(tmp>eps) {
               if(tmp>err) err=tmp;
               diff=fii-fjj;
               sign=1.0;
               if (diff<0) {
                  diff=-diff;
                  sign=-sign;
               }
               d=diff+sqrt(diff*diff+4*fij*fij);
               t=2*sign*fij/d;
               c=1/sqrt(1+t*t);
               s=c*t;
/* rotate u */
               for (k=0; k<nr; k++) {
                  d=-s*u[k+i*nr]+c*u[k+j*nr];
                  u[k+i*nr]=c*u[k+i*nr]+s*u[k+j*nr];
                  u[k+j*nr]=d;
               }
/* rotate f with on index in i or j */
               for (k=0; k<n; k++) {
                  if ( (k!=i) || (k!=j) ) {
                     ki=k+i*(i+1)/2;
                     kj=k+j*(j+1)/2;
                     if (k>i) ki=i+k*(k+1)/2;
                     if (k>j) kj=j+k*(k+1)/2;
                     d=c*f[kj]-s*f[ki];
                     f[ki]=s*f[kj]+c*f[ki];
                     f[kj]=d;
                  }
               }
               f[ii]=c*c*fii+s*s*fjj+2*c*s*fij;
               f[jj]=s*s*fii+c*c*fjj-2*c*s*fij;
               f[ij]=0;
            }
         }
      }
   }
}
