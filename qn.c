/****************************************************************************/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "expopt.h"
#define GUESS
#define NEWLS
#define MAX_ITER 50
#define BFGS      1
#define DPF       2
#define GS        1.61803398875
typedef struct {
   double   energy;
   double  *x;
   double  *g;
   double  *h;
} point;
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void quasi_newton(double x[], int n, double ethr, double gthr, FILE *f) {
   double *B,*Bd;
   point   P[MAX_ITER+1];
   double *t;
   double *s;
   double *h;
#ifndef NEWLS
   double  gamma,hBig;
#endif
   double  ener[99],step[99];
   double  e_0,e_1,e_est,e_p,e_m;
   double  d,sg,ss,lambda;
   double  x1,x2,g1,g2;
   double  gd,dBd;
   double  gnorm;
   double  trust;
   double  p;
   int     iter,k0;
   int     np;
   int     i,j,k;
   fprintf(f,"Entering quasi_newton\n");
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   trust=0.50;
/*--------------------------------------------------------------------------*/
/* Allocate.                                                                */
/*--------------------------------------------------------------------------*/
   for(k=0; k<=MAX_ITER; k++) {
      P[k].x=(double *)malloc(n*sizeof(double));
      P[k].g=(double *)malloc(n*sizeof(double));
      P[k].h=(double *)malloc(n*sizeof(double));
   }
   B=(double *)malloc(n*n*sizeof(double));
   Bd=(double *)malloc(n*sizeof(double));
   t=(double *)malloc(n*sizeof(double));
   s=(double *)malloc(n*sizeof(double));
   h=(double *)malloc(n*sizeof(double));
/*--------------------------------------------------------------------------*/
/* Initialize.                                                              */
/*--------------------------------------------------------------------------*/
   for(i=0; i<n; i++) P[0].x[i]=x[i];
   for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
         if(i==j) B[i+n*j]=SYS.B0;
         else     B[i+n*j]=0.0;
      }
   }
   for(i=0; i<n; i++) h[i]=0.10;
/*--------------------------------------------------------------------------*/
/* Iterate.                                                                 */
/*--------------------------------------------------------------------------*/
   fprintf(f,"Iter            E(0)            E(est)          E(1)           dE          gnorm        l.s.   \n");
   fprintf(f,"Iter       -------------   -------------   -------------  -----------   ----------   ----------\n");
   iter=0;
   e_0=energy(P[0].x);
   while(iter<MAX_ITER) {
      for(i=0; i<n; i++) t[i]=P[iter].x[i];
      P[iter].energy=e_0;
      fprintf(f,"e_0=%.8f\n",e_0);
      for(j=0; j<n; j++) {
         d=sqrt(1.0e-5/B[j+n*j]);
         t[j]=t[j]+d;
         e_p=energy(t);
         t[j]=t[j]-2.0*d;
         e_m=energy(t);
         t[j]=t[j]+d;
         fprintf(f,"e_p[%i]=%.8f\n",j,e_p);
         fprintf(f,"e_m[%i]=%.8f\n",j,e_m);
         P[iter].g[j]=(e_p-e_m)/(2.0*d);
         P[iter].h[j]=(e_p-2.0*e_0+e_m)/(d*d);
      }
      fprintf(f,"Gradient:");
      for(i=0; i<n; i++) fprintf(f,"%10.6f ",P[iter].g[i]); fprintf(f,"\n");
      fprintf(f,"Hessian: ");
      for(i=0; i<n; i++) fprintf(f,"%10.6f ",P[iter].h[i]); fprintf(f,"\n");
      gnorm=0.0;
      for(i=0; i<n; i++) gnorm=gnorm+P[iter].g[i]*P[iter].g[i];
      gnorm=sqrt(gnorm);
      for(i=0; i<n; i++) h[i]=sqrt(2.0e-5/P[iter].h[i]);

      k0=iter-10;
      if(k0<0) k0=0;
      for(k=k0; k<=iter; k++) {
         if(k==k0) {
#ifdef GUESS
            guessian(B,P[k0].h,n,f);
#else
            hBig=0.0;
            for(i=0; i<n; i++) if(hBig<P[k0].h[i]) hBig=P[k0].h[i];
            for(i=0; i<n; i++) {
               p=P[k0].h[i];
               if(p<1.0e-3*hBig) p=1.0e-3*hBig;
               for(j=0; j<n; j++) B[i+n*j]=0.0;
               B[i+n*i]=p;
            }
#endif
         } else {
            for(i=0; i<n; i++) {
               p=0.0;
               for(j=0; j<n; j++) p=p+B[i+n*j]*(P[k].x[j]-P[k-1].x[j]);
               Bd[i]=p;
            }
	    printf("g:  "); for(i=0; i<n; i++) printf("%10.6f ",(P[k].g[i]-P[k-1].g[i])); printf("\n");
	    printf("d:  "); for(i=0; i<n; i++) printf("%10.6f ",(P[k].x[i]-P[k-1].x[i])); printf("\n");
	    printf("Bd: "); for(i=0; i<n; i++) printf("%10.6f ",Bd[i]); printf("\n");
            gd=0.0;
            dBd=0.0;
            for(i=0; i<n; i++) {
               gd=gd+(P[k].g[i]-P[k-1].g[i])*(P[k].x[i]-P[k-1].x[i]);
               dBd=dBd+(P[k].x[i]-P[k-1].x[i])*Bd[i];
            }
	    printf("gd: %10.6f\n",gd);
	    printf("dBd: %10.6f\n",dBd);
            for(i=0; i<n; i++) {
               for(j=0; j<n; j++) {
                  B[i+n*j]=B[i+n*j]
                          +(P[k].g[i]-P[k-1].g[i])*(P[k].g[j]-P[k-1].g[j])/gd
                          -Bd[i]*Bd[j]/dBd;
               }
            }
         }
         fprintf(f,"Hessian:\n");
         for(i=0; i<n; i++) {
            for(j=0; j<n; j++) fprintf(f,"%10.6f ",B[i+n*j]); fprintf(f,"\n");
         }
      }
      if(1) check_hessian(B,n);
      gauss(B,s,P[iter].g,n);
      for(i=0; i<n; i++) s[i]=-s[i];
      fprintf(f,"Step:    ");
      for(i=0; i<n; i++) fprintf(f,"%10.6f ",s[i]); fprintf(f,"\n");

#ifdef NEWLS
      sg=0.0; for(i=0; i<n; i++) sg=sg+P[iter].g[i]*s[i];
      ss=0.0; for(i=0; i<n; i++) ss=ss+s[i]*s[i];
      e_est=e_0+0.5*sg;
      fprintf(f,"e_est=%.8f\n",e_est);

      step[0]=0.0;
      ener[0]=e_0;
      np=1;

      if(trust/sqrt(ss)>1) lambda=1.0;
      else                 lambda=trust/sqrt(ss);
      fprintf(f,"lambda=%.6f\n",lambda);

      if(0) {
         step[1]=lambda/(1.0+GS);
      } else if(1) {
         step[1]=lambda/(1.0+GS);
         if(step[1]<0.002/sqrt(ss)) step[1]=0.002/sqrt(ss);
         if(step[1]>0.050/sqrt(ss)) step[1]=0.050/sqrt(ss);
      } else {
         step[1]=0.1;
      }
      for(i=0; i<n; i++) P[iter+1].x[i]=P[iter].x[i]+step[1]*s[i];
      ener[1]=energy(P[iter+1].x);
      np=2;
      
      if(ener[1]<ener[0]) {
         step[2]=(1.0+GS)*step[1]-GS*step[0];
         for(i=0; i<n; i++) P[iter+1].x[i]=P[iter].x[i]+step[2]*s[i];
         ener[2]=energy(P[iter+1].x);
         np=3;
      } else {
         fprintf(f,"qn: FUNNNNNNNY, wrong way!?\n");
         p      =step[0];
         step[0]=step[1];
         step[1]=p;
         p      =ener[0];
         ener[0]=ener[1];
         ener[1]=p;
         step[2]=2.0*step[1]-1.0*step[0];
         for(i=0; i<n; i++) P[iter+1].x[i]=P[iter].x[i]+step[2]*s[i];
         ener[2]=energy(P[iter+1].x);
         np=3;
      }

      for(i=0; i<np; i++) {
         fprintf(f,"step[%2i]=%.8f, ",i,step[i]);
         fprintf(f,"ener[%2i]=%.8f",i,ener[i]);
         fprintf(f,"\n");
      }

      k=find_min(ener,np);
      while(k==np-1) {
         step[np]=(1.0+GS)*step[np-1]-GS*step[np-2];
         for(i=0; i<n; i++) P[iter+1].x[i]=P[iter].x[i]+step[np]*s[i];
         ener[np]=energy(P[iter+1].x);
         fprintf(f,"step[%2i]=%.8f, ",np,step[np]);
         fprintf(f,"ener[%2i]=%.8f",np,ener[np]);
         fprintf(f,"\n");
         np++;
         k=find_min(ener,np);
      }

      k=find_min(ener,np);
      if(k==np-1) {
         lambda=step[k];
      } else {
         x1=0.5*(step[k]+step[k-1]);
         x2=0.5*(step[k+1]+step[k]);
         g1=(ener[k]-ener[k-1])/(step[k]-step[k-1]);
         g2=(ener[k+1]-ener[k])/(step[k+1]-step[k]);
         lambda=(x1*g2-x2*g1)/(g2-g1);
      }
      fprintf(f,"lambda=%.6f\n",lambda);

      for(i=0; i<n; i++) P[iter+1].x[i]=P[iter].x[i]+lambda*s[i];
      e_1=energy(P[iter+1].x);
      fprintf(f,"e_1=%.8f\n",e_1);
      if(e_1>ener[k]+1.0e-10) {
         lambda=step[k];
         for(i=0; i<n; i++) P[iter+1].x[i]=P[iter].x[i]+lambda*s[i];
         e_1=ener[k];
      }
      
#else
      for(i=0; i<n; i++) P[iter+1].x[i]=P[iter].x[i]+s[i];
      e_1=energy(P[iter+1].x);
      fprintf(f,"e_1=%.8f\n",e_1);

      sg=0.0;
      for(i=0; i<n; i++) sg=sg+P[iter].g[i]*s[i];
      e_est=e_0+0.5*sg;
      fprintf(f,"e_est=%.8f\n",e_est);

      if( (e_0-e_est) > 1.10*(e_0-e_1) ) {
         if(0) {
            gamma=3.0*(e_est-e_1)/sg;
            fprintf(f,"gamma=%.8f\n",gamma);
            lambda=2.0/(1.0+sqrt(1.0+4.0*gamma));
            fprintf(f,"lambda=%.8f\n",lambda);
         } else {
            gamma=2.0*(e_est-e_1)/sg;
            fprintf(f,"gamma=%.8f\n",gamma);
            lambda=1.0/(1.0+gamma);
            fprintf(f,"lambda=%.8f\n",lambda);
         }
         for(i=0; i<n; i++) P[iter+1].x[i]=P[iter].x[i]+lambda*s[i];
         e_1=energy(P[iter+1].x);
         fprintf(f,"e_1=%.8f\n",e_1);
      }
#endif
      for(i=0; i<n; i++) x[i]=P[iter+1].x[i];

      fprintf(f,"Iter %2i: %15.8f %15.8f %15.8f %12.8f %12.8f %12.8f\n",iter,e_0,e_est,e_1,e_1-e_0,gnorm,lambda);
      if( ((e_0-e_1)<ethr) && (gnorm<gthr) ) {
         fprintf(f,"Converged\n");
         break;
      };
      e_0=e_1;
      iter++;
   }
   e_0=e_1;
/*--------------------------------------------------------------------------*/
/* Deallocate vertices.                                                     */
/*--------------------------------------------------------------------------*/
   for(k=0; k<=MAX_ITER; k++) {
      free(P[k].x);
      free(P[k].g);
      free(P[k].h);
   }
   free(B);
   free(Bd);
   free(t);
   free(s);
   free(h);
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   fprintf(f,"Exiting quasi_newton\n");
}
