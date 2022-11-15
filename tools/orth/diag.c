/****************************************************************************/
/*                                                                          */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "onbas.h"
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void diag(basis *Y) {
   basis   X;
   double *t_prim;
   double *v_prim;
   double *h_prim;
   double *c_prim;
   double *t_cont;
   double *v_cont;
   double *h_cont;
   double *c_cont;
   double *e_cont;
   double  ft,fv,x1,x2;
   double  alpha,beta;
   double  t;
   int     mxprim;
   int     mxcont;
   int     nprim;
   int     ncont;
   int     mxlqn;
   int     lqn;
   int     ind;
   int     i,j,k,m;
#ifdef TRACE
   fprintf(stderr,"+++ Entering diag\n");
#endif
   X=*Y;
   mxlqn=X.maxlqn;
/*--------------------------------------------------------------------------*/
/* Determine dimensions and allocate                                        */
/*--------------------------------------------------------------------------*/
   mxprim=0;
   mxcont=0;
   for(lqn=0; lqn<=mxlqn; lqn++) {
      if(mxprim<X.shell[lqn].nprim) mxprim=X.shell[lqn].nprim;
      if(mxcont<X.shell[lqn].ncont) mxcont=X.shell[lqn].ncont;
   }
#ifdef DEBUG
   fprintf(stderr,"mxprim=%i\n",mxprim);
   fprintf(stderr,"mxcont=%i\n",mxcont);
#endif
   t_prim=(double *)malloc(mxprim*mxprim*sizeof(double));
   v_prim=(double *)malloc(mxprim*mxprim*sizeof(double));
   h_prim=(double *)malloc(mxprim*mxprim*sizeof(double));
   c_prim=(double *)malloc(mxprim*mxcont*sizeof(double));
   t_cont=(double *)malloc((mxcont*(mxcont+1))/2*sizeof(double));
   v_cont=(double *)malloc((mxcont*(mxcont+1))/2*sizeof(double));
   h_cont=(double *)malloc((mxcont*(mxcont+1))/2*sizeof(double));
   c_cont=(double *)malloc(mxcont*mxcont*sizeof(double));
   e_cont=(double *)malloc(mxcont*sizeof(double));
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   fv=-2.0*X.Z/sqrt(PI);
   for(lqn=0; lqn<=mxlqn; lqn++) {
#ifdef DEBUG
      fprintf(stderr,"Diagonalizing lqn=%i\n",lqn);
#endif
      nprim=X.shell[lqn].nprim;
      ncont=X.shell[lqn].ncont;
/*
 Construct primitive integrals
*/
      for(i=0; i<nprim; i++) {
         alpha=X.shell[lqn].z[i];
         for(j=0; j<nprim; j++) {
            beta=X.shell[lqn].z[j];
            ind=i+nprim*j;
            ft=0.25*(2.0*(double)(lqn)+3.0);
            x1=0.25*(2.0*(double)(lqn)+7.0);
            x2=0.50*(2.0*(double)(lqn)+5.0);
            t_prim[ind]=ft*pow(4.0*alpha*beta,x1)/pow(alpha+beta,x2);
            x1=0.25*(2.0*(double)(lqn)+3.0);
            x2=(double)(lqn)+1.0;
            v_prim[ind]=fv*pow(4.0*alpha*beta,x1)/pow(alpha+beta,x2);
            h_prim[ind]=t_prim[ind]+v_prim[ind];
         }
      }
#ifdef DEBUG
      fprintf(stderr,"t_prim\n");
      for(i=0; i<nprim; i++) {
         for(j=0; j<nprim; j++) {
            ind=i+nprim*j;
            fprintf(stderr,"%12.6f ",t_prim[ind]);
         }
         fprintf(stderr,"\n");
      }
      fprintf(stderr,"v_prim\n");
      for(i=0; i<nprim; i++) {
         for(j=0; j<nprim; j++) {
            ind=i+nprim*j;
            fprintf(stderr,"%12.6f ",v_prim[ind]);
         }
         fprintf(stderr,"\n");
      }
      fprintf(stderr,"h_prim\n");
      for(i=0; i<nprim; i++) {
         for(j=0; j<nprim; j++) {
            ind=i+nprim*j;
            fprintf(stderr,"%12.6f ",h_prim[ind]);
         }
         fprintf(stderr,"\n");
      }
#endif
      fv=fv*(2.0*(double)(lqn)+2.0)/(2.0*(double)(lqn)+3.0);
/*
 Construct contracted integrals
*/
      for(k=0; k<ncont; k++) for(m=0; m<=k; m++) {
         t=0.0;
         for(i=0; i<nprim; i++) for(j=0; j<nprim; j++) t=t+X.shell[lqn].c[i][k]*t_prim[i+nprim*j]*X.shell[lqn].c[j][m];
         t_cont[k*(k+1)/2+m]=t;
      }
      for(k=0; k<ncont; k++) for(m=0; m<=k; m++) {
         t=0.0;
         for(i=0; i<nprim; i++) for(j=0; j<nprim; j++) t=t+X.shell[lqn].c[i][k]*v_prim[i+nprim*j]*X.shell[lqn].c[j][m];
         v_cont[k*(k+1)/2+m]=t;
      }
      for(k=0; k<ncont; k++) for(m=0; m<=k; m++) {
         t=0.0;
         for(i=0; i<nprim; i++) for(j=0; j<nprim; j++) t=t+X.shell[lqn].c[i][k]*h_prim[i+nprim*j]*X.shell[lqn].c[j][m];
         h_cont[k*(k+1)/2+m]=t;
      }
      for(k=0; k<ncont; k++) { for(m=0; m<ncont; m++) c_cont[m+ncont*k]=0.0; c_cont[k+ncont*k]=1.0; }
#ifdef DEBUG
      fprintf(stderr,"t_cont\n");
      ind=0; for(k=0; k<ncont; k++) { for(m=0; m<=k; m++) { fprintf(stderr,"%12.6f ",t_cont[ind]); ind++; } fprintf(stderr,"\n"); }
      fprintf(stderr,"v_cont\n");
      ind=0; for(k=0; k<ncont; k++) { for(m=0; m<=k; m++) { fprintf(stderr,"%12.6f ",v_cont[ind]); ind++; } fprintf(stderr,"\n"); }
      fprintf(stderr,"h_cont\n");
      ind=0; for(k=0; k<ncont; k++) { for(m=0; m<=k; m++) { fprintf(stderr,"%12.6f ",h_cont[ind]); ind++; } fprintf(stderr,"\n"); }
      fprintf(stderr,"c_cont\n");
      ind=0; for(k=0; k<ncont; k++) { for(m=0; m<ncont; m++) { fprintf(stderr,"%12.6f ",c_cont[ind]); ind++; } fprintf(stderr,"\n"); }
#endif
/*
 Diagonalize one electron hamiltonian
*/
      jacobi(h_cont,c_cont,ncont,ncont);
      for(k=0; k<ncont; k++) e_cont[k]=h_cont[k*(k+1)/2+k];
#ifdef DEBUG
      fprintf(stderr,"h_cont after diagonalization\n");
      ind=0; for(k=0; k<ncont; k++) { for(m=0; m<=k; m++) { fprintf(stderr,"%12.6f ",h_cont[ind]); ind++; } fprintf(stderr,"\n"); }
      fprintf(stderr,"c_cont after diagonalization\n");
      ind=0; for(k=0; k<ncont; k++) { for(m=0; m<ncont; m++) { fprintf(stderr,"%12.6f ",c_cont[ind]); ind++; } fprintf(stderr,"\n"); }
#endif
      jacord(e_cont,c_cont,ncont,ncont);
#ifdef DEBUG
      fprintf(stderr,"eigenvalues and eigenvectors\n");
      for(k=0; k<ncont; k++) fprintf(stderr,"%12.6f ",e_cont[k]); fprintf(stderr,"\n");
      for(k=0; k<ncont; k++) fprintf(stderr,"------------ "); fprintf(stderr,"\n");
      ind=0; for(k=0; k<ncont; k++) { for(m=0; m<ncont; m++) { fprintf(stderr,"%12.6f ",c_cont[ind]); ind++; } fprintf(stderr,"\n"); }
#endif
/*
 Transform basis set
*/
      for(i=0; i<nprim; i++) {
         for(k=0; k<ncont; k++) {
            t=0.0;
            for(m=0; m<ncont; m++) {
               t=t+X.shell[lqn].c[i][m]*c_cont[m+k*ncont];
            }
            c_prim[i+k*nprim]=t;
         }
      }
#ifdef DEBUG
      fprintf(stderr,"c_prim\n");
      for(i=0; i<nprim; i++) { for(k=0; k<ncont; k++) fprintf(stderr,"%12.6f ",c_prim[i+k*nprim]); fprintf(stderr,"\n"); }
#endif
      for(i=0; i<nprim; i++) {
         for(k=0; k<ncont; k++) {
            X.shell[lqn].c[i][k]=c_prim[i+k*nprim];
         }
      }
   }
   *Y=X;
/*--------------------------------------------------------------------------*/
/* Deallocate                                                               */
/*--------------------------------------------------------------------------*/
   free(e_cont);
   free(c_cont);
   free(h_cont);
   free(v_cont);
   free(t_cont);
   free(c_prim);
   free(h_prim);
   free(v_prim);
   free(t_prim);
/*--------------------------------------------------------------------------*/
/* Done                                                                     */
/*--------------------------------------------------------------------------*/
#ifdef TRACE
   fprintf(stderr,"+++ Exiting diag\n");
#endif
}
