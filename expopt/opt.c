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
void optimize(int l, int ind[], int n, FILE *f) {
   double *x;
   char    shell[]={'s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u','v','w','x','y','z'};
   int     i,k;
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   fprintf(f,"***\n");
   fprintf(f,"*** Optimizing shell %i\n",l);
   fprintf(f,"***\n");
   if(SYS.generator==GEN_HOT) {
      fprintf(f,"   Hot tempered generator\n");
      fprintf(f,"   Directions: "); for(k=0; k<n; k++) fprintf(f,"%i ",ind[k]); fprintf(f,"\n");
   } else if(SYS.generator==GEN_HOTPO) {
      fprintf(f,"   Hot tempered generator with polynomial extension\n");
      fprintf(f,"   Directions: "); for(k=0; k<n; k++) fprintf(f,"%i ",ind[k]); fprintf(f,"\n");
   } else if(SYS.generator==GEN_HARM) {
      fprintf(f,"   Harmonic generator\n");
      fprintf(f,"   Directions: "); for(k=0; k<n; k++) fprintf(f,"%i ",ind[k]); fprintf(f,"\n");
   } else if(SYS.generator==GEN_POLY) {
      fprintf(f,"   Polynomial generator\n");
      fprintf(f,"   Directions: "); for(k=0; k<n; k++) fprintf(f,"%i ",ind[k]); fprintf(f,"\n");
   } else if(SYS.generator==GEN_COS) {
      fprintf(f,"   Cosine generator\n");
      fprintf(f,"   Directions: "); for(k=0; k<n; k++) fprintf(f,"%i ",ind[k]); fprintf(f,"\n");
   } else if(SYS.generator==GEN_UNIT) {
      fprintf(f,"   Unit generator\n");
      fprintf(f,"   Directions: "); for(k=0; k<n; k++) fprintf(f,"%i ",ind[k]); fprintf(f,"\n");
   } else if(SYS.generator==GEN_EXP) {
      fprintf(f,"   Exponential generator, alpha=%.6f\n",SYS.alpha);
      fprintf(f,"   Directions: "); for(k=0; k<n; k++) fprintf(f,"%i ",ind[k]); fprintf(f,"\n");
   } else if(SYS.generator==GEN_SPLIT) {
      fprintf(f,"   Split valence generator\n");
      fprintf(f,"   Directions: "); for(k=0; k<n; k++) fprintf(f,"%i ",ind[k]); fprintf(f,"\n");
   } else if(SYS.generator==GEN_INPUT) {
      fprintf(f,"   Input directions\n");
   } else {
      fprintf(stderr,"Error: Internal error in routine optimize\n");
      fprintf(stderr,"       Illegal generator type!\n");
      exit(1);
   }
   if(fabs(SYS.pen_a)>0.10e-6) {
      fprintf(f,"   Penalty function %.8f %.8f %.8f\n",SYS.pen_a,SYS.pen_b,SYS.pen_c);
   }
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   SYS.energy=1.0e12;
   SYS.active=l;
   x=(double *)malloc(n*sizeof(double));
   for(k=0; k<n; k++) x[k]=0.0;
   if(SYS.index!=NULL) free(SYS.index);
   SYS.index=(int *)malloc(n*sizeof(int));
   for(k=0; k<n; k++) SYS.index[k]=ind[k];
   SYS.n_dir=n;

   if(n==1) {
      powell(x,n,SYS.ethr,f);
      fprintf(f,"After powell: ");
      for(k=0; k<n; k++) {fprintf(f,"%12.6f ",x[k]);} fprintf(f,"\n");
   } else if(SYS.method==OPT_AUTO) {
      fprintf(stderr,"Error: Internal error in routine optimize\n");
      fprintf(stderr,"       Method 'Auto' not yet implemented!\n");
      exit(1);
   } else if(SYS.method==OPT_AMOEBA) {
      amoeba(x,n,SYS.ethr,f);
      fprintf(f,"After amoeba: ");
      for(k=0; k<n; k++) {fprintf(f,"%12.6f ",x[k]);} fprintf(f,"\n");
   } else if(SYS.method==OPT_POWELL) {
      powell(x,n,SYS.ethr,f);
      fprintf(f,"After powell: ");
      for(k=0; k<n; k++) {fprintf(f,"%12.6f ",x[k]);} fprintf(f,"\n");
   } else if(SYS.method==OPT_QN) {
      quasi_newton(x,n,SYS.ethr,SYS.gthr,f);
      fprintf(f,"After quasi newton: ");
      for(k=0; k<n; k++) {fprintf(f,"%12.6f ",x[k]);} fprintf(f,"\n");
   } else if(SYS.method==OPT_QN2) {
/*      qn2(x,n,SYS.ethr,SYS.gthr,f); */
      fprintf(f,"qn2 is not yet implemented: ");
      for(k=0; k<n; k++) {fprintf(f,"%12.6f ",x[k]);} fprintf(f,"\n");
      exit(1);
   } else {
      fprintf(stderr,"Error: Internal error in routine optimize\n");
      fprintf(stderr,"       Illegal method type!\n");
      exit(1);
   }

   for(i=0; i<SYS.basis[l].n; i++) SYS.basis[l].z[i]=SYS.basis[l].z[i]*make_scale(i,x);
   fprintf(f,"Optimized exponents:\n");
   if(SYS.generator==GEN_HOT) {
      fprintf(f,"%s %i%c: HT/%i, E=%.8f\n","%",SYS.basis[l].n,shell[l],n,SYS.energy);
   } else if(SYS.generator==GEN_HOTPO) {
      fprintf(f,"%s %i%c: HP/%i, E=%.8f\n","%",SYS.basis[l].n,shell[l],n,SYS.energy);
   } else if(SYS.generator==GEN_HARM) {
      fprintf(f,"%s %i%c: HG/%i, E=%.8f\n","%",SYS.basis[l].n,shell[l],n,SYS.energy);
   } else if(SYS.generator==GEN_POLY) {
      fprintf(f,"%s %i%c: PG/%i, E=%.8f\n","%",SYS.basis[l].n,shell[l],n,SYS.energy);
   } else if(SYS.generator==GEN_COS) {
      fprintf(f,"%s %i%c: CG/%i, E=%.8f\n","%",SYS.basis[l].n,shell[l],n,SYS.energy);
   } else if(SYS.generator==GEN_UNIT) {
      fprintf(f,"%s %i%c: UG/%i, E=%.8f\n","%",SYS.basis[l].n,shell[l],n,SYS.energy);
   } else if(SYS.generator==GEN_EXP) {
      fprintf(f,"%s %i%c: EG/%i(%.6f), E=%.8f\n","%",SYS.basis[l].n,shell[l],n,SYS.alpha,SYS.energy);
   } else if(SYS.generator==GEN_SPLIT) {
      fprintf(f,"%s %i%c: SV/%i, E=%.8f\n","%",SYS.basis[l].n,shell[l],n,SYS.energy);
   } else if(SYS.generator==GEN_INPUT) {
      fprintf(f,"%s %i%c: IN/%i, E=%.8f\n","%",SYS.basis[l].n,shell[l],n,SYS.energy);
   } else {
      fprintf(stderr,"Error: Internal error in routine optimize\n");
      fprintf(stderr,"       Illegal generator type!\n");
      exit(1);
   }
   fprintf(f," %2i",SYS.basis[l].Z);
   fprintf(f," %2i",SYS.basis[l].l);
   fprintf(f," %2i ",SYS.basis[l].n);
   for(i=0; i<SYS.basis[l].n; i++) {
      if( (i!=0) && (i%5==0) ) fprintf(f,"\n          ");
      fprintf(f," %.8f",SYS.basis[l].z[i]);
   }
   fprintf(f,"\n");

   free(SYS.index); SYS.index=NULL;
   free(x); x=NULL;
   SYS.active=-1;
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
}
