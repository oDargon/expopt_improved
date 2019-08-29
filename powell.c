/****************************************************************************/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "expopt.h"
#define MAX_ITER  100
#define MAX_DIR    32
#define GS        1.61803398875
typedef struct {
   double  *x;
   double   h;
} DIR;
typedef struct {
   double  *x;
   double   e;
} POINT;
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void powell(double z[], int n, double thr, FILE *f) {
   DIR     D[MAX_DIR+1];
   POINT   P[MAX_DIR+1];
   double  x[MAX_DIR];
   double  s[MAX_DIR];
   double  M[MAX_DIR][MAX_DIR];
   double  e0,e1;
   double  p,q;
   int     iter;
   int     i,j,k,m;
   fprintf(f,"Entering powell\n");
/*--------------------------------------------------------------------------*/
/* Initialize.                                                              */
/*--------------------------------------------------------------------------*/
   fprintf(f,"thr=%.10f\n",thr);
   if(n>MAX_DIR) {
      fprintf(stderr,"Maximum number of directione exceeded in powell!\n");
      exit(1);
   }
/* prime directions */
   for(k=0; k<n+1; k++) {
      D[k].x=(double *)malloc(n*sizeof(double));
      for(i=0; i<n; i++) D[k].x[i]=0.0;
      D[k].x[k]=1.0;
      printf("***\n");
      printf("*** Fix initial step length in powell!!!\n");
      printf("***\n");
      if(0) D[k].h=0.150;
      else  D[k].h=0.004;
   }
/* prime points */
   for(k=0; k<MAX_DIR+1; k++) {
      P[k].x=(double *)malloc(n*sizeof(double));
      for(i=0; i<n; i++) P[0].x[i]=0.0;
   }
/* starting point */
   for(i=0; i<n; i++) x[i]=z[i];
   for(i=0; i<n; i++) P[0].x[i]=x[i];
/*--------------------------------------------------------------------------*/
/* Iterate.                                                                 */
/*--------------------------------------------------------------------------*/
   iter=0;
   e0=P[0].e=energy(P[0].x);
   fprintf(f,"Iter           Energy          dE              x               x               x\n");
   fprintf(f,"Iter        ------------    ------------    ------------    ------------    ------------\n");
   fprintf(f,"Iter     %15.8f\n",P[0].e);
   while(iter<MAX_ITER) {
      fprintf(f,"\n");
      fprintf(f,"...\n");
      fprintf(f,"... Iteration %i\n",iter);
      fprintf(f,"...\n");
      for(k=0; k<n+1; k++) { /* loop over all directions */
         fprintf(f,"Line search in direction %2i (%.8f)",k,D[k].h);
         for(i=0; i<n; i++) fprintf(f," %12.6f",D[k].x[i]); fprintf(f,"\n");
         if(fabs(D[k].h)<1.0e-4) {
            if(D[k].h<0.0) D[k].h=-1.0e-4;
            else           D[k].h= 1.0e-4;
         }
         for(i=0; i<n; i++) s[i]=D[k].x[i];
         if(n==1) e1=line_search(x,s,n,&D[k].h,1.0e-6,f);
         else     e1=line_search(x,s,n,&D[k].h,1.0e-6,f); /* quick fix */
         if(n==1) break;
         D[k].h=0.3*D[k].h;
         P[k+1].e=e1;
         for(i=0; i<n; i++) P[k+1].x[i]=x[i];
         if(k==n-1) { /* Assemble the composite direction */
            fprintf(f,"Assembling composite direction!\n");
            for(i=0; i<n; i++) D[n].x[i]=P[n].x[i]-P[0].x[i];
            D[n].h=0.1;
            fprintf(f,"Direction %2i:",n);
            for(i=0; i<n; i++) fprintf(f," %12.6f",D[n].x[i]);
            fprintf(f,"\n");
            p=0.0;
            for(i=0; i<n; i++) p=p+D[n].x[i]*D[n].x[i];
            p=1.0/sqrt(p);
            for(i=0; i<n; i++) D[n].x[i]=p*D[n].x[i];
            D[n].h=0.3/p;
            fprintf(f,"Direction %2i:",n);
            for(i=0; i<n; i++) fprintf(f," %12.6f",D[n].x[i]);
            fprintf(f,"\n");
         }
      }
      if(n==1) break;
      if(1) for(k=0; k<n+2; k++) { /* Debug print points */
         if(k==0) fprintf(f,"Point %2i: %15.8f %15.8f",k,P[k].e,0.0);
         else     fprintf(f,"Point %2i: %15.8f %15.8f",k,P[k].e,P[k].e-P[k-1].e);
         for(i=0; i<n; i++) fprintf(f," %12.6f",P[k].x[i]);
         fprintf(f,"\n");
      }
      p=0.0; q=0.0;
      for(i=0; i<n; i++) { p=p + (P[n+1].x[i]-P[n  ].x[i]) * (P[n  ].x[i]-P[0].x[i]); fprintf(f,"%.6f,",p); } fprintf(f,"\n");
      for(i=0; i<n; i++) { q=q + (P[n  ].x[i]-P[0  ].x[i]) * (P[n  ].x[i]-P[0].x[i]); fprintf(f,"%.6f,",q); } fprintf(f,"\n");
      fprintf(f,"p/q=%.10f\n",p/q);
      if(p/q>0.5) { /* Replace the direction with the largest decrease */
         m=0; for(k=0; k<n; k++) if(P[k+1].e-P[k].e<P[m+1].e-P[m].e) m=k;
         fprintf(f,"Big=%i\n",m);
         for(i=0; i<n; i++) D[m].x[i]=D[n].x[i];
         D[m].h=D[n].h;
         for(k=0; k<n; k++) {
            fprintf(f,"Dir %2i: %12.6f",k,D[k].h);
            for(i=0; i<n; i++) fprintf(f," %12.6f",D[k].x[i]);
            fprintf(f,"\n");
         }
         if(1) { /* Check that we have linear independence */
            for(i=0; i<n; i++) for(j=0; j<n; j++) M[i][j]=D[i].x[j];
            p=check_lin_dep(M,n,f);
            fprintf(f,"Smallest eigenvalue %.6f\n",p);
            if(p<1.0e-4) {
               fprintf(f,"... Replacing directions!\n");
               for(k=0; k<n+1; k++) {
                  for(i=0; i<n; i++) D[k].x[i]=M[i][k];
                  D[k].h=0.05;
               }
            }
         }
      }
      fprintf(f,"Iter %2i: %15.8f %15.8f %15.8f %15.8f %15.8f\n",iter,e1,e1-e0,0.0,0.0,0.0);
      P[0].e=P[n].e;
      for(i=0; i<n; i++) P[0].x[i]=P[n].x[i];
      if(fabs(e0-e1)<thr) break;
      e0=e1;
      iter++;
   } /* End iterations */
   for(i=0; i<n; i++) z[i]=x[i];
/*--------------------------------------------------------------------------*/
/* Deallocate vertices.                                                     */
/*--------------------------------------------------------------------------*/
   for(k=0; k<MAX_DIR+1; k++) {
      free(P[k].x);
   }
   for(k=0; k<n+1; k++) {
      free(D[k].x);
   }
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   fprintf(f,"Exiting powell\n");
}
