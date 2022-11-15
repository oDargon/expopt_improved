/****************************************************************************/
/*                                                                          */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "onbas.h"
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void rdbas(FILE *f, basis *Y) {
   basis  X;
   char   line[256];
   char  *tok;
   int    lqn;
   int    i,k;
#ifdef TRACE
   fprintf(stderr,"+++ Entering rdbas\n");
#endif
   if(fgets(line,sizeof(line),f)==NULL) { fprintf(stderr,"Premature end of input!\n"); exit(1); }
   tok=strtok(line," \t\n"); X.Z=atof(tok);
   tok=strtok(NULL," \t\n"); X.maxlqn=atoi(tok);
   for(lqn=0; lqn<=X.maxlqn; lqn++) {
#ifdef DEBUG
      fprintf(stderr,"Reading lqn=%i\n",lqn);
#endif
      if(fgets(line,sizeof(line),f)==NULL) { fprintf(stderr,"Premature end of input!\n"); exit(1); }
      tok=strtok(line," \t\n"); X.shell[lqn].nprim=atoi(tok);
      tok=strtok(NULL," \t\n"); X.shell[lqn].ncont=atoi(tok);
      tok=NULL;
      for(k=0; k<X.shell[lqn].nprim; k++) {
         if(tok==NULL) {
            if(fgets(line,sizeof(line),f)==NULL) { fprintf(stderr,"Premature end of input!\n"); exit(1); }
            tok=strtok(line," \t\n");
         }
         X.shell[lqn].z[k]=atof(tok);
         tok=strtok(NULL," \t\n");
      }
      for(k=0; k<X.shell[lqn].nprim; k++) {
         if(fgets(line,sizeof(line),f)==NULL) { fprintf(stderr,"Premature end of input!\n"); exit(1); }
         tok=strtok(line," \t\n");
         for(i=0; i<X.shell[lqn].ncont; i++) {
            X.shell[lqn].c[k][i]=atof(tok);
            tok=strtok(NULL," \t\n");
         }
      }
   if(fgets(line,sizeof(line),f)==NULL) { fprintf(stderr,"Premature end of input!\n"); exit(1); }
   tok=strtok(line," \t\n"); X.En=atof(tok);
   if(X.En>0){
         if(fgets(line,sizeof(line),f)==NULL) { fprintf(stderr,"Premature end of input!\n"); exit(1); }
         }
   }
#ifdef DEBUG
   fprintf(stderr,"Z=%.0f, maxlqn=%i\n",X.Z,X.maxlqn);
   for(lqn=0; lqn<=X.maxlqn; lqn++) {
      fprintf(stderr,"nprim[%i]=%i, ncont[%i]=%i\n",lqn,X.shell[lqn].nprim,lqn,X.shell[lqn].ncont);
   }
#endif
   *Y=X;
#ifdef TRACE
   fprintf(stderr,"+++ Exiting rdbas\n");
#endif
}
