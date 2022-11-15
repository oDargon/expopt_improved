/****************************************************************************/
/*                                                                          */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "onbas.h"
#define DEBUG
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void gram(basis *Y) {
   basis  X;
   double zi,zj;
   double ci,cj;
   double p,q,t;
   int    lqn;
   int    ic,jc;
   int    ip,jp;
   double THR=0.000001;
#ifdef TRACE
   fprintf(stderr,"+++ Entering gram\n");
#endif
   X=*Y;
#ifdef DEBUG
   fprintf(stderr,"Z=%.0f, maxlqn=%i\n",X.Z,X.maxlqn);
   for(lqn=0; lqn<=X.maxlqn; lqn++) {
      fprintf(stderr,"nprim[%i]=%i, ncont[%i]=%i\n",lqn,X.shell[lqn].nprim,lqn,X.shell[lqn].ncont);
   }
#endif
   p=0.75;
   for(lqn=0; lqn<=X.maxlqn; lqn++) {
#ifdef DEBUG
      {
         double S[MXCONT][MXCONT];
         fprintf(stderr,"Orthogonalizing lqn=%i\n",lqn);
         for(ic=0; ic<X.shell[lqn].ncont; ic++) {
            for(jc=0; jc<X.shell[lqn].ncont; jc++) {
               t=0.0;
               for(ip=0; ip<X.shell[lqn].nprim; ip++) {
                  zi=X.shell[lqn].z[ip];
                  ci=X.shell[lqn].c[ip][ic];
                  for(jp=0; jp<X.shell[lqn].nprim; jp++) {
                     zj=X.shell[lqn].z[jp];
                     cj=X.shell[lqn].c[jp][jc];
                     q=4.0*zi*zj/((zi+zj)*(zi+zj));
                     t=t+ci*cj*pow(q,p);
                  }
               }
               S[ic][jc]=t;
            }
         }
         for(ic=0; ic<X.shell[lqn].ncont; ic++) {
            for(jc=0; jc<X.shell[lqn].ncont; jc++) fprintf(stderr,"%12.6f ",S[ic][jc]);
            fprintf(stderr,"\n");
         }
         
         for(ic=0;ic<X.shell[lqn].ncont; ic++) {
         for(jc=0; jc<X.shell[lqn].ncont; jc++){
           if(ic==jc && (fabs(S[ic][jc]-1.0)>THR)) { fprintf (stderr,"Kaboom Z=%f shell=%d element=%d,%d %18.12f\n",X.Z,lqn,ic,ic,S[ic][jc]); exit(1);}
           if(ic!=jc && (fabs(S[ic][jc])>THR)) { fprintf (stderr,"Kaboom Z=%f shell=%d element=%d,%d %18.12f\n", X.Z,lqn,ic,jc,S[ic][jc]); exit(1); }
         }
         }
      }
#endif
      for(ic=0; ic<X.shell[lqn].ncont; ic++) {
         /* Normalize */
         t=0.0;
         for(ip=0; ip<X.shell[lqn].nprim; ip++) {
            zi=X.shell[lqn].z[ip];
            ci=X.shell[lqn].c[ip][ic];
            for(jp=0; jp<X.shell[lqn].nprim; jp++) {
               zj=X.shell[lqn].z[jp];
               cj=X.shell[lqn].c[jp][ic];
               q=4.0*zi*zj/((zi+zj)*(zi+zj));
               t=t+ci*cj*pow(q,p);
            }
         }
         t=1.0/sqrt(t);
         for(ip=0; ip<X.shell[lqn].nprim; ip++) X.shell[lqn].c[ip][ic]=t*X.shell[lqn].c[ip][ic];
         /* Orthogonalize */
         for(jc=ic+1; jc<X.shell[lqn].ncont; jc++) {
            t=0.0;
            for(ip=0; ip<X.shell[lqn].nprim; ip++) {
               zi=X.shell[lqn].z[ip];
               ci=X.shell[lqn].c[ip][ic];
               for(jp=0; jp<X.shell[lqn].nprim; jp++) {
                  zj=X.shell[lqn].z[jp];
                  cj=X.shell[lqn].c[jp][jc];
                  q=4.0*zi*zj/((zi+zj)*(zi+zj));
                  t=t+ci*cj*pow(q,p);
               }
            }
            for(ip=0; ip<X.shell[lqn].nprim; ip++) X.shell[lqn].c[ip][jc]=X.shell[lqn].c[ip][jc]-t*X.shell[lqn].c[ip][ic];
         }
      }
      p=p+0.5;
   }
   *Y=X;
#ifdef TRACE
   fprintf(stderr,"+++ Exiting gram\n");
#endif
}
