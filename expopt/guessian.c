/****************************************************************************/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "expopt.h"
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void guessian(double B[], double diag[], int nd, FILE *f) {
   double *H;
   double *V;
   double *s;
   double  dBig;
   double  t;
   int     np,l;
   int     i,j,p,q;
   fprintf(f,"Entering guessian\n");
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   l=SYS.active;
   np=SYS.basis[l].n;
/*--------------------------------------------------------------------------*/
/* Allocate.                                                                */
/*--------------------------------------------------------------------------*/
   H=(double *)malloc(np*np*sizeof(double));
   V=(double *)malloc(nd*np*sizeof(double));
   s=(double *)malloc(nd*sizeof(double));
/*--------------------------------------------------------------------------*/
/* Build guessed primitive hessian.                                         */
/*--------------------------------------------------------------------------*/
   for(i=0; i<np; i++) {
      for(j=i; j<np; j++) {
	 if(i==j)             { H[i+np*j]= 1.0; H[j+np*i]= 1.0; }
	 else if(abs(i-j)==1) { H[i+np*j]=-0.7; H[j+np*i]=-0.7; }
	 else if(abs(i-j)==2) { H[i+np*j]= 0.3; H[j+np*i]= 0.3; }
	 else if(abs(i-j)==3) { H[i+np*j]=-0.1; H[j+np*i]=-0.1; }
	 else                 { H[i+np*j]= 0.0; H[j+np*i]= 0.0; }
      }
   }
   fprintf(f,"Guessed primitive hessian\n");
   for(i=0; i<np; i++) {
      for(j=0; j<np; j++) {fprintf(f,"%7.4f ",H[i+np*j]);} fprintf(f,"\n");
   }
/*--------------------------------------------------------------------------*/
/* Build transformation vectors.                                            */
/*--------------------------------------------------------------------------*/
   for(p=0; p<nd; p++) {
      for(q=0; q<nd; q++) s[q]=0.0;
      s[p]=1.0;
      for(i=0; i<np; i++) {
         t=make_scale(i,s);
         V[i+np*p]=log(t);
      }
   }
   fprintf(f,"Transformation matrix\n");
   for(i=0; i<np; i++) {
      for(p=0; p<nd; p++) {fprintf(f,"%7.4f ",V[i+np*p]);} fprintf(f,"\n");
   }
/*--------------------------------------------------------------------------*/
/* Transform hessian. Done as O(n4), can be done as O(n3).                  */
/*--------------------------------------------------------------------------*/
   for(p=0; p<nd; p++) {
      for(q=0; q<nd; q++) {
         t=0.0;
         for(i=0; i<np; i++) {
            for(j=0; j<np; j++) {
               t=t+V[i+np*p]*V[j+np*q]*H[i+np*j];
            }
         }
         B[p+nd*q]=t;
      }
   }
   fprintf(f,"Transformed hessian\n");
   for(p=0; p<nd; p++) {
      for(q=0; q<nd; q++) {fprintf(f,"%7.4f ",B[p+nd*q]);} fprintf(f,"\n");
   }
/*--------------------------------------------------------------------------*/
/* Scale transformed hessian.                                               */
/*--------------------------------------------------------------------------*/
   dBig=diag[0]; for(p=0; p<nd; p++) if(diag[p]>dBig) dBig=diag[p];
   if(dBig<1.0e-4) dBig=1.0e-4;
   for(p=0; p<nd; p++) {
      t=diag[p];
      if(t<5.0e-3*dBig) t=5.0e-3*dBig;
      s[p]=sqrt(t/B[p+nd*p]);
   }
   for(p=0; p<nd; p++) {
      for(q=0; q<nd; q++) {
         B[p+nd*q]=B[p+nd*q]*s[p]*s[q];
      }
   }
   fprintf(f,"Scaled transformed hessian\n");
   for(p=0; p<nd; p++) {
      for(q=0; q<nd; q++) {fprintf(f,"%7.4f ",B[p+nd*q]);} fprintf(f,"\n");
   }
/*--------------------------------------------------------------------------*/
/* Deallocate.                                                              */
/*--------------------------------------------------------------------------*/
   free(H);
   free(V);
   free(s);
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   fprintf(f,"Exiting guessian\n");
}
