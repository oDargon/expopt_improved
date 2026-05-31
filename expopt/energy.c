/****************************************************************************/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "expopt.h"
double energy(double *s) {
   FILE  *f;
   char   line[128];
   int    k,l,m,n;
   int    max_lqn;
   double e,z,p;
   e=0.0;
/*--------------------------------------------------------------------------*/
/* Print basis set into to file.                                            */
/*--------------------------------------------------------------------------*/
   if((f=fopen("BASIS.DAT","w"))==NULL) {
      fprintf(stderr,"Could not open file BASIS.DAT\n");
      exit(1);
   }
   max_lqn=-1;
   for(l=0; l<=MAX_LQN; l++) {
      if(SYS.basis[l].n>0) max_lqn=l;
   }
   if(SYS.printthis==-1) {
      fprintf(f,"%6.2f %i\n",SYS.charge,max_lqn);
      for(l=0; l<=max_lqn; l++) {
         n=SYS.basis[l].n;
         if(n==0) {
            fprintf(f," 0 0\n");
         } else {
            int frozen = SYS.contract_frozen && (l != SYS.active);
            int nc     = frozen ? SYS.basis[l].m : n;
            fprintf(f," %i %i\n",n,nc);
            for(k=0; k<n; k++) {
               z=SYS.basis[l].z[k];
               if(SYS.active==l) {
                  z=z*make_scale(k,s);
               }
               fprintf(f,"%.8f\n",z);
            }
            if(frozen) {
               for(k=0; k<nc; k++) {
                  for(m=0; m<n; m++) fprintf(f,"%.8f ",SYS.basis[l].c[k*n+m]);
                  fprintf(f,"\n");
               }
            } else {
               for(k=0; k<n; k++) {
                  for(m=0; m<n; m++) {
                     if(k==m) fprintf(f,"1. ");
                     else     fprintf(f,"0. ");
                  }
                  fprintf(f,"\n");
               }
            }
         }
      }
   } else if(0) {
      l=SYS.printthis;
      n=SYS.basis[l].n;
      for(k=0; k<n; k++) {
         z=SYS.basis[l].z[k];
         if(SYS.active==l) {
            z=z*make_scale(k,s);
         }
         fprintf(f,"%.8f\n",z);
      }
   } else {
      l=SYS.printthis;
      n=SYS.basis[l].n;
      if(SYS.printraw) {
         for(k=0; k<n; k++) {
            z=SYS.basis[l].z[k];
            if(SYS.active==l) {
               z=z*make_scale(k,s);
            }
            fprintf(f," %.8f",z);
         }
         fprintf(f,"\n");
      } else {
         if(n==0) {
            fprintf(f," 0 0\n");
         } else {
            fprintf(f," %i %i\n",n,n);
            for(k=0; k<n; k++) {
               z=SYS.basis[l].z[k];
               if(SYS.active==l) {
                  z=z*make_scale(k,s);
               }
               fprintf(f,"%.8f\n",z);
            }
            for(k=0; k<n; k++) {
               for(m=0; m<n; m++) {
                  if(k==m) fprintf(f,"1. ");
                  else     fprintf(f,"0. ");
               }
               fprintf(f,"\n");
            }
         }
      }
   }
   fclose(f);
/*--------------------------------------------------------------------------*/
/* Execute energy module.                                                   */
/*--------------------------------------------------------------------------*/
   fflush(stdout);
   if(system("./energy.ksh")!=0) {
      fprintf(stderr,"Error executing 'energy'\n");
      exit(1);
   }
/*--------------------------------------------------------------------------*/
/* Read energy from file.                                                   */
/*--------------------------------------------------------------------------*/
   if((f=fopen("ENERGY.DAT","r"))==NULL) {
      fprintf(stderr,"Could not open file 'ENERGY.DAT'\n");
      exit(1);
   }
   if(fgets(line,sizeof(line),f)==NULL) {
      fprintf(stderr,"Premature end of 'ENERGY.DAT'\n");
      exit(1);
   }
   e=atof(line);
   fclose(f);
/*--------------------------------------------------------------------------*/
/* Add penalty function to energy                                           */
/*--------------------------------------------------------------------------*/
   if( (SYS.active>-1) ) {
      p=0.0;
      l=SYS.active;
      n=SYS.basis[l].n;
      for(k=0; k<n; k++) {
         z=SYS.basis[l].z[k];
         z=z*make_scale(k,s);
         z=log(z);
         p=p+SYS.pen_a*exp(SYS.pen_b*(z-SYS.pen_c));
      }
      e=e+p;
   }
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   if( (SYS.active>-1) && (e<SYS.energy) ) {
      l=SYS.active;
      n=SYS.basis[l].n;
      if(1) {
         printf("Energy %.8f\n",e);
         printf(" %2i ",SYS.basis[l].Z);
         printf("%2i ",SYS.basis[l].l);
         printf("%2i  ",SYS.basis[l].n);
         for(k=0; k<n; k++) {
            if( (k!=0) && (k%5==0) ) printf("\n           ");
            z=SYS.basis[l].z[k];
            z=z*make_scale(k,s);
            printf("%.8f ",z);
         }
         printf("\n");
      }
      SYS.energy=e;
   }
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   SYS.fnc_eval++;
   return e;
}
