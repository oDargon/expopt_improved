/****************************************************************************/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "expopt.h"
typedef struct {
   double *x;
} dir;
double make_scale(int K, double *s) {
   dir    *P,Q;
   double  p,q,f,x;
   double  alpha;
   double  pi;
   int     n_prim,n_dir;
   int     i,l,k,m;
/*--------------------------------------------------------------------------*/
/* Allocate directions.                                                     */
/*--------------------------------------------------------------------------*/
   l=SYS.active;
   n_prim=SYS.basis[l].n;
   n_dir=SYS.n_dir;
   P=(dir *)malloc(n_dir*sizeof(dir));
   for(k=0; k<n_dir; k++) P[k].x=(double *)malloc(n_prim*sizeof(double));
   Q.x=(double *)malloc(n_prim*sizeof(double));
/*--------------------------------------------------------------------------*/
/* Generate basis directions.                                               */
/*--------------------------------------------------------------------------*/
   if(SYS.generator==GEN_HOT) {
      for(k=0; k<n_dir; k++) {
         if(SYS.index[k]==0) {
            for(i=0; i<n_prim; i++) P[k].x[i]=1.0;
         } else if(SYS.index[k]==1) {
            for(i=0; i<n_prim; i++) P[k].x[i]=(double)(i)-0.5*(double)(n_prim-1);
         } else {
            q=1.0-(double)(SYS.index[k]);
            for(i=0; i<n_prim; i++) P[k].x[i]=pow((double)(i+1),q);
         }
      }
   } else if(SYS.generator==GEN_HOTPO) {
      for(k=0; k<n_dir; k++) {
         if(SYS.index[k]==0) {
            fprintf(stderr,"DEBUG: direction x^0\n");
            for(i=0; i<n_prim; i++) P[k].x[i]=1.0;
         } else if(SYS.index[k]==1) {
            fprintf(stderr,"DEBUG: direction x^1\n");
            for(i=0; i<n_prim; i++) P[k].x[i]=(double)(i)-0.5*(double)(n_prim-1);
         } else if(SYS.index[k]==2) {
            fprintf(stderr,"DEBUG: direction x^{-1}\n");
            for(i=0; i<n_prim; i++) P[k].x[i]=1.0/(double)(i+1);
         } else {
            q=(double)(SYS.index[k])-1.0;
            fprintf(stderr,"DEBUG: direction x^{%.0f}\n",q);
            for(i=0; i<n_prim; i++) P[k].x[i]=pow((double)(i+1),q);
         }
      }
   } else if(SYS.generator==GEN_SPLIT) {
      for(k=0; k<n_dir; k++) {
         if(SYS.index[k]==0) {
            for(i=0; i<n_prim; i++) P[k].x[i]=1.0;
         } else if(SYS.index[k]==1) {
            for(i=0; i<n_prim; i++) P[k].x[i]=(double)(i)-0.5*(double)(n_prim-1);
         } else if(SYS.index[k]==2) {
            for(i=0; i<n_prim; i++) P[k].x[i]=1.0/(double)(i+1);
         } else {
            for(i=0; i<n_prim; i++) P[k].x[i]=0.0;
            i=n_prim+2-SYS.index[k];
            P[k].x[i]=1.0;
         }
      }
   } else if(SYS.generator==GEN_EXP) {
      alpha=SYS.alpha;
      for(k=0; k<n_dir; k++) {
         if(SYS.index[k]==0) {
            for(i=0; i<n_prim; i++) P[k].x[i]=1.0;
         } else if(SYS.index[k]==1) {
            for(i=0; i<n_prim; i++) P[k].x[i]=(double)(i)-0.5*(double)(n_prim-1);
         } else if(SYS.index[k]==2) {
            for(i=0; i<n_prim; i++) P[k].x[i]=exp(-alpha*(double)(i));
         } else if(SYS.index[k]==3) {
            for(i=0; i<n_prim; i++) P[k].x[i]=(double)(i)*exp(-alpha*(double)(i));
         } else if(SYS.index[k]==4) {
            for(i=0; i<n_prim; i++) P[k].x[i]=exp(-alpha*(double)(n_prim-1-i));
         } else if(SYS.index[k]==5) {
            for(i=0; i<n_prim; i++) P[k].x[i]=(double)(n_prim-1-i)*exp(-alpha*(double)(n_prim-1-i));
         } else {
            fprintf(stderr,"Error: illegal direction specified for exponential direction\n");
            exit(1);
         }
      }
   } else if(SYS.generator==GEN_HARM) {
      fprintf(stderr,"Error: internal error in make_scale\n");
      fprintf(stderr,"       Harmonic generator not implemented!\n");
      exit(1);
   } else if(SYS.generator==GEN_POLY) { /* Chebychev polynomials */
      for(k=0; k<n_dir; k++) {
         for(i=0; i<n_prim; i++) {
            x=-1.0+2.0*(double)(i)/(double)(n_prim-1);
            x=0.999999*x;
            P[k].x[i]=cos((double)(SYS.index[k])*acos(x));
         }
      }
   } else if(SYS.generator==GEN_COS) {
      pi=4.0*atan(1.0);
      for(k=0; k<n_dir; k++) {
         for(i=0; i<n_prim; i++) {
            P[k].x[i]=cos(pi*(double)(i)*(double)(SYS.index[k])/(double)(n_prim-1));
         }
      }
   } else if(SYS.generator==GEN_UNIT) {
      for(k=0; k<n_dir; k++) {
         for(i=0; i<n_prim; i++) P[k].x[i]=0.0;
         P[k].x[SYS.index[k]]=1.0;
      }
   } else if(SYS.generator==GEN_INPUT) {
      for(k=0; k<n_dir; k++) {
         for(i=0; i<n_prim; i++) {
            P[k].x[i]=SYS.directions[i+k*n_prim];
         }
      }
   } else {
      fprintf(stderr,"Error: internal error in make_scale\n");
      fprintf(stderr,"       Illegal generator!\n");
      exit(1);
   }
#ifdef DEBUG
   printf("Directions\n");
   for(k=0; k<n_dir; k++) {
      for(i=0; i<n_prim; i++) printf("%12.6f ",P[k].x[i]);
      printf("\n");
   }
#endif
/*--------------------------------------------------------------------------*/
/* Orthonormalize directions.                                               */
/*--------------------------------------------------------------------------*/
   for(k=0; k<n_dir; k++) {
      for(m=0; m<k; m++) {
         p=0.0;
         for(i=0; i<n_prim; i++) p=p+P[k].x[i]*P[m].x[i];
         for(i=0; i<n_prim; i++) P[k].x[i]=P[k].x[i]-p*P[m].x[i];
      }
      p=0.0;
      for(i=0; i<n_prim; i++) p=p+P[k].x[i]*P[k].x[i];
      p=sqrt(1.0/p);
      for(i=0; i<n_prim; i++) P[k].x[i]=p*P[k].x[i];
   }
#ifdef DEBUG
   printf("Directions\n");
   for(k=0; k<n_dir; k++) {
      for(i=0; i<n_prim; i++) printf("%12.6f ",P[k].x[i]);
      printf("\n");
   }
   printf("Overlaps:\n");
   for(k=0; k<n_dir; k++) {
      for(m=0; m<=k; m++) {
         p=0.0;
         for(i=0; i<n_prim; i++) p=p+P[k].x[i]*P[m].x[i];
         printf("%12.6f ",p);
      }
      printf("\n");
   }
#endif
/*--------------------------------------------------------------------------*/
/* Assemble step                                                            */
/*--------------------------------------------------------------------------*/
   for(i=0; i<n_prim; i++) {
      p=0.0;
      for(k=0; k<n_dir; k++) {
         p=p+s[k]*P[k].x[i];
      }
      Q.x[i]=p;
   }
/*--------------------------------------------------------------------------*/
/* Compute scale factor                                                     */
/*--------------------------------------------------------------------------*/
   f=exp(Q.x[K]);
#ifdef DEBUG
   printf("Scale factor=%.6f\n",f);
#endif
/*--------------------------------------------------------------------------*/
/* Deallocate directions.                                                   */
/*--------------------------------------------------------------------------*/
   for(k=0; k<n_dir; k++) free(P[k].x);
   free(P);
   free(Q.x);
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   return f;
}
