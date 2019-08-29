/****************************************************************************/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "expopt.h"
typedef struct {
   double   energy;
   double  *x;
} point;
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void amoeba(double x[], int n, double thr, FILE *f) {
   point *P,Q;
   double e_hi,e_lo;
   double step;
   double t;
   int    iter,hi,hi2,lo;
   int    i,j,k,m;
   fprintf(f,"Entering amoeba\n");
/*--------------------------------------------------------------------------*/
/* Allocate vertices.                                                       */
/*--------------------------------------------------------------------------*/
   P=(point *)malloc((n+1)*sizeof(point));
   for(k=0; k<=n; k++) P[k].x=(double *)malloc(n*sizeof(double));
   Q.x=(double *)malloc(n*sizeof(double));
/*--------------------------------------------------------------------------*/
/* Construct vertices with unit distances.                                  */
/*--------------------------------------------------------------------------*/
   for(k=0; k<=n; k++) {
      for(i=0; i<n; i++) P[k].x[i]=0.0;
   }
   P[1].x[0]=1.0;
   for(k=2; k<=n; k++) {
      for(i=0; i<n; i++) {
         t=0.0;
         for(j=0; j<k; j++) t=t+P[j].x[i];
         P[k].x[i]=t/(double)(k);
      }
      t=0.0;
      for(i=0; i<n; i++) t=t+P[k].x[i]*P[k].x[i];
      t=sqrt(1.0-t);
      P[k].x[k-1]=t;
   }
   fprintf(f,"Vertices\n");
   for(k=0; k<=n; k++) {
      for(i=0; i<n; i++) fprintf(f,"%12.6f ",P[k].x[i]);
      fprintf(f,"\n");
   }
   fprintf(f,"Distances\n");
   for(k=1; k<=n; k++) {
      for(m=0; m<k; m++) {
         t=0.0;
         for(i=0; i<n; i++) t=t+(P[k].x[i]-P[m].x[i])*(P[k].x[i]-P[m].x[i]);
         t=sqrt(t);
         fprintf(f,"%12.6f ",t);
      }
      printf("\n");
   }
/*--------------------------------------------------------------------------*/
/* Transform vertices.                                                      */
/*--------------------------------------------------------------------------*/
   step=SYS.A0;
   fprintf(f,"Initial step in amoeba %.6f\n",step);
   for(k=0; k<=n; k++) {
      for(i=0; i<n; i++) {
         P[k].x[i]=step*P[k].x[i]+x[i];
      }
   }
   fprintf(f,"Vertices\n");
   for(k=0; k<=n; k++) {
      for(i=0; i<n; i++) fprintf(f,"%12.6f ",P[k].x[i]);
      fprintf(f,"\n");
   }
   fprintf(f,"Distances\n");
   for(k=1; k<=n; k++) {
      for(m=0; m<k; m++) {
         t=0.0;
         for(i=0; i<n; i++) t=t+(P[k].x[i]-P[m].x[i])*(P[k].x[i]-P[m].x[i]);
         t=sqrt(t);
         fprintf(f,"%12.6f ",t);
      }
      printf("\n");
   }
/*--------------------------------------------------------------------------*/
/* Prime vertices.                                                          */
/*--------------------------------------------------------------------------*/
   for(k=0; k<=n; k++) {
      P[k].energy=energy(P[k].x);
   }
   fprintf(f,"Primed vertices\n");
   for(k=0; k<=n; k++) {
      for(i=0; i<n; i++) fprintf(f,"%12.6f ",P[k].x[i]);
      fprintf(f,"%16.10f\n",P[k].energy);
   }
/*--------------------------------------------------------------------------*/
/* Iterate.                                                                 */
/*--------------------------------------------------------------------------*/
   for(iter=0; iter<250; iter++) {
      printf(">>> Iteration %i\n",iter+1);

      if(1) {
         fprintf(f,"Vertices\n");
         for(k=0; k<=n; k++) {
            if(1) for(i=0; i<n; i++) fprintf(f,"%12.6f ",P[k].x[i]);
            fprintf(f,"%16.10f\n",P[k].energy);
         }
      }
      if(0) {
         fprintf(f,"Distances\n");
         for(k=1; k<=n; k++) {
            for(m=0; m<k; m++) {
               t=0.0;
               for(i=0; i<n; i++) t=t+(P[k].x[i]-P[m].x[i])*(P[k].x[i]-P[m].x[i]);
               t=sqrt(t);
               fprintf(f,"%12.6f ",t);
            }
            printf("\n");
         }
      }

      hi=0; 
      for(k=1; k<=n; k++) if(P[k].energy>P[hi].energy) hi=k;
      printf("Vertex %2i is high with energy        %.10f\n",hi,P[hi].energy);

      if(hi==0) hi2=1;
      else      hi2=0;
      for(k=0; k<=n; k++) if( (k!=hi) && (P[k].energy>P[hi2].energy) ) hi2=k;
      printf("Vertex %2i is second high with energy %.10f\n",hi2,P[hi2].energy);

      lo=0;
      for(k=1; k<=n; k++) if(P[k].energy<P[lo].energy) lo=k;
      printf("Vertex %2i is low with energy         %.10f\n",lo,P[lo].energy);

      for(i=0; i<n; i++) {
         t=0.0;
         for(k=0; k<=n; k++) {
            if(k!=hi) {
               t=t+2.0*P[k].x[i]-1.0*P[hi].x[i];
            }
         }
         Q.x[i]=t/(double)(n);
      }
      Q.energy=energy(Q.x);
      fprintf(f,"New vertex with energy               %.10f\n",Q.energy);
      for(i=0; i<n; i++) fprintf(f,"%12.6f ",Q.x[i]);
      fprintf(f,"\n");

      e_hi=P[hi].energy;
      e_lo=P[lo].energy;
      printf("fabs(e_hi-e_lo)=%.10f\n",fabs(e_hi-e_lo));
      if(Q.energy<P[hi2].energy) {
         printf("Amoeba: reflection!\n");
         for(i=0; i<n; i++) P[hi].x[i]=Q.x[i];
         P[hi].energy=Q.energy;
      } else if(fabs(e_hi-e_lo)>thr) {
         printf("Amoeba: shrinking!\n");
         if(Q.energy<P[hi].energy) {
            for(i=0; i<n; i++) P[hi].x[i]=Q.x[i];
            P[hi].energy=Q.energy;
         }
         t=sqrt(0.1);
         for(k=0; k<=n; k++) {
            if(k!=lo) {
               for(i=0; i<n; i++) P[k].x[i]=t*P[k].x[i]+(1.0-t)*P[lo].x[i];
               P[k].energy=energy(P[k].x);
            }
         }
      } else {
         break;
      }
 
   }
   if(Q.energy<P[lo].energy) {
      printf("Q is low\n");
      for(i=0; i<n; i++) x[i]=Q.x[i];
   } else {
      printf("P[%i] is low\n",lo);
      for(i=0; i<n; i++) x[i]=P[lo].x[i];
   }
/*--------------------------------------------------------------------------*/
/* Deallocate vertices.                                                     */
/*--------------------------------------------------------------------------*/
   for(k=0; k<=n; k++) free(P[k].x);
   free(P);
   free(Q.x);
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   fprintf(f,"Exiting amoeba\n");
}
