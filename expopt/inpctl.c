/****************************************************************************/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "expopt.h"
#define NEWLINE 1
#define MUSTHAVE 2
/*--------------------------------------------------------------------------*/
#define GETLINE \
if(fgets(line,sizeof(line),iu)==NULL) { \
   fprintf(stderr,"Unexpected end of input file\n"); \
   exit(1); \
} \
line_no++;
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void input_control(FILE *iu, FILE *ou) {
   char    line[4096],*tok;
   double  e;
   int     ind[128];
   int     line_no;
   int     Z,l,n,m;
   int     k;
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   line_no=0;
   GETLINE
   tok=strtok(line," \n");
   if(strcmp(tok,"#ExpOpt-2.0")!=0) {
      fprintf(stderr,"Wrong file type!\n");
      exit(1);
   }
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   while(fgets(line,sizeof(line),iu)!=NULL) {
      line_no++;
      if(SYS.verbose) fprintf(ou,"+++ %s",line);
      if(line[0]=='#') continue;
      tok=strtok(line," \n");
/*
   Define basis set
*/
      if(strcmp(tok,"Define")==0) {
         GETLINE
         tok=strtok(line," \n"); Z=atoi(tok);
         tok=strtok(NULL," \n"); l=atoi(tok);
         tok=strtok(NULL," \n"); n=atoi(tok);
         tok=strtok(NULL," \n"); m=atoi(tok);
         if(l>MAX_LQN) {
            fprintf(stderr,"Max. l-q.n. exceeded\n");
            exit(1);
         }
         SYS.basis[l].Z=Z;
         SYS.basis[l].l=l;
         SYS.basis[l].n=n;
         SYS.basis[l].m=m;
         if(SYS.basis[l].z!=NULL) {
            fprintf(ou,"Freeing shell %i\n",l);
            free(SYS.basis[l].z);
         }
         if(SYS.basis[l].c!=NULL) free(SYS.basis[l].c);
         SYS.basis[l].z=(double *)malloc(n*sizeof(double));
         SYS.basis[l].c=(double *)malloc(m*n*sizeof(double));
         for(k=0; k<n; k++) {
            tok=strtok(NULL," \n");
            while(tok==NULL) {
               GETLINE
               tok=strtok(line," \n");
            }
            SYS.basis[l].z[k]=atof(tok);
         }
         for(k=0; k<m*n; k++) {
            tok=strtok(NULL," \n");
            while(tok==NULL) {
               GETLINE
               tok=strtok(line," \n");
            }
            SYS.basis[l].c[k]=atof(tok);
         }
/*
   Print current basis set
*/
      } else if(strcmp(tok,"Print")==0) {
         for(l=0; l<=MAX_LQN; l++) {
            if(SYS.basis[l].n>0) {
               fprintf(ou,"--- Shell %i\n",l);
               fprintf(ou,"    %2i ",SYS.basis[l].Z);
               fprintf(ou,"%2i ",SYS.basis[l].l);
               fprintf(ou,"%2i ",SYS.basis[l].n);
               fprintf(ou,"%2i ",SYS.basis[l].m);
               for(k=0; k<SYS.basis[l].n; k++) {
                  if( (k!=0) && ((k%5)==0) ) fprintf(ou,"\n             ");
                  fprintf(ou,"%.6f ",SYS.basis[l].z[k]);
               }
               fprintf(ou,"\n");
               for(m=0; m<SYS.basis[l].m; m++) {
                  fprintf(ou,"    ");
                  for(k=0; k<SYS.basis[l].n; k++) fprintf(ou,"%.6f ",SYS.basis[l].c[m*SYS.basis[l].n+k]);
                  fprintf(ou,"\n");
               }
            }
         }
/*
   Compute energy
*/
      } else if(strcmp(tok,"Energy")==0) {
         e=energy(NULL);
         fprintf(ou,"Energy = %.8f\n",e);
/*
   Set charge
*/
      } else if(strcmp(tok,"Charge")==0) {
         tok=strtok(NULL," \n"); SYS.charge=atof(tok);
/*
   Optimize shell
*/
      } else if(strcmp(tok,"Optimize")==0) {
         tok=strtok(NULL," \n"); l=atoi(tok);
         n=0; 
         while((tok=strtok(NULL," \n"))!=NULL) {
            ind[n]=atoi(tok);
            n++;
         }
         if(n<1) {
            fprintf(stderr,"Too few directions specified on line %i\n",line_no);
            exit(1);
         }
         optimize(l,ind,n,ou);
/*
   Optimize shell with input directions
*/
      } else if(strcmp(tok,"OptimizeWithDirections")==0) {
         SYS.generator=GEN_INPUT;
         tok=strtok(NULL," \n"); l=atoi(tok);
         tok=strtok(NULL," \n"); n=atoi(tok);
         tok=NULL;
         for(k=0; k<n*SYS.basis[l].n; k++) {
            if(tok==NULL) {
               GETLINE
               tok=strtok(line," \n");
            }
            SYS.directions[k]=atof(tok);
            tok=strtok(NULL," \n");
         }
         /*
         for(k=0; k<n*SYS.basis[l].n; k++) {
            if( (k!=0) && (k%SYS.basis[l].n==0) ) printf("\n");
            printf("%6.3f ",SYS.directions[k]);
         }
         printf("\n");
         */
         optimize(l,ind,n,ou);
/*
   Print search directions
*/
      } else if(strcmp(tok,"PrintDirections")==0) {
         fprintf(ou,"!!! PrintDirections\n");
         exit(1);
/*
   Line search shell
*/
      } else if(strcmp(tok,"LineSearch")==0) {
         fprintf(ou,"!!! LineSearch\n");
         exit(1);
/*
   Compute hessian
*/
      } else if(strcmp(tok,"Hessian")==0) {
         fprintf(ou,"!!! Hessian\n");
         exit(1);
/*
   Select method
*/
      } else if(strcmp(tok,"Method")==0) {
         tok=strtok(NULL," \n");
         if(strcmp(tok,"Auto")==0) {
            SYS.method=OPT_AUTO;
         } else if(strcmp(tok,"Amoeba")==0) {
            SYS.method=OPT_AMOEBA;
         } else if(strcmp(tok,"Powell")==0) {
            SYS.method=OPT_POWELL;
         } else if(strcmp(tok,"QN")==0) {
            SYS.method=OPT_QN;
/*         } else if(strcmp(tok,"QN2")==0) {
            SYS.method=OPT_QN2; */
         } else {
            fprintf(stderr,"Illegal method '%s'\n",tok);
            exit(1);
         }
/*
   Select generator
*/
      } else if(strcmp(tok,"Generator")==0) {
         tok=strtok(NULL," \n");
         if(strcmp(tok,"HotTempered")==0) {
            SYS.generator=GEN_HOT;
         } else if(strcmp(tok,"Polynomial")==0) {
            SYS.generator=GEN_POLY;
         } else if(strcmp(tok,"Cosine")==0) {
            SYS.generator=GEN_COS;
         } else if(strcmp(tok,"Unit")==0) {
            SYS.generator=GEN_UNIT;
         } else if(strcmp(tok,"Exponential")==0) {
            tok=strtok(NULL," \n"); SYS.alpha=atof(tok);
            SYS.generator=GEN_EXP;
         } else if(strcmp(tok,"SplitValence")==0) {
            SYS.generator=GEN_SPLIT;
         } else if(strcmp(tok,"HotPoly")==0) {
            SYS.generator=GEN_HOTPO;
         } else {
            fprintf(stderr,"Illegal generator type '%s'\n",tok);
            exit(1);
         }
/*
   Select to output only selected shell.
*/
      } else if(strcmp(tok,"SelectShell")==0) {
         tok=strtok(NULL," \n"); SYS.printthis=atoi(tok);
/*
   Select to output only selected shell.
*/
      } else if(strcmp(tok,"PrintRaw")==0) {
         SYS.printraw=1;
/*
   Contract frozen shells during optimization
*/
      } else if(strcmp(tok,"ContractFrozen")==0) {
         SYS.contract_frozen=1;
/*
   Skip line search in QN, use trust-region step only
*/
      } else if(strcmp(tok,"NoLineSearch")==0) {
         SYS.no_linesearch=1;
/*
   Set the penalty function
*/
      } else if(strcmp(tok,"Penalty")==0) {
         tok=strtok(NULL," \n"); SYS.pen_a=atof(tok);
         tok=strtok(NULL," \n"); SYS.pen_b=atof(tok);
         tok=strtok(NULL," \n"); SYS.pen_c=atof(tok);
/*
   Change Z for all shells
*/
      } else if(strcmp(tok,"ChangeZ")==0) {
         tok=strtok(NULL," \n"); Z=atoi(tok);
         for(k=0; k<=MAX_LQN; k++) SYS.basis[k].Z=Z;
/*
   Set initial step length
*/
      } else if(strcmp(tok,"Step")==0) {
         double d;
         tok=strtok(NULL," \n"); d=atof(tok);
         SYS.B0=1.0e-5/(d*d);
         SYS.A0=d;
/*
   Sort the primitives
*/
      } else if(strcmp(tok,"Sort")==0) {
         fprintf(ou,"!!! Sort\n");
         exit(1);
/*
   Reorder the primitives
*/
      } else if(strcmp(tok,"Reorder")==0) {
         fprintf(ou,"!!! Reorder\n");
         exit(1);
/*
   Fit the primitives
*/
      } else if(strcmp(tok,"Fit")==0) {
         fprintf(ou,"!!! Fit\n");
         exit(1);
/*
   Set energy threshold
*/
      } else if(strcmp(tok,"Threshold")==0) {
         tok=strtok(NULL," \n");
         SYS.ethr=atof(tok);
         if(SYS.ethr<1.0e-10) SYS.ethr=1.0e-10;
/*
   Set gradient threshold
*/
      } else if(strcmp(tok,"gThreshold")==0) {
         tok=strtok(NULL," \n");
         SYS.gthr=atof(tok);
         if(SYS.gthr<1.0e-5) SYS.gthr=1.0e-5;
/*
   Extend shell
*/
      } else if(strcmp(tok,"Extend")==0) {
         tok=strtok(NULL," \n"); l=atof(tok);
         tok=strtok(NULL," \n"); n=atof(tok);
         extend(l,n,ou);
/*
   Error: keyword not understood
*/
      } else {
         fprintf(stderr,"Keyword '%s' ???\n",tok);
         exit(1);
      }
   }
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
}
