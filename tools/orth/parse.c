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
void parse(int argc, char *argv[]) {
   char parm[128];
   int  k,m;
#ifdef TRACE
   fprintf(stderr,"+++ Entering parse\n");
#endif
   ORTHO=ORTHO_GRAM;
   DIAG=0;
   k=1;
   while(k<argc) {
      strcpy(parm,argv[k]);
      if(parm[0]!='-') {
         syntax(stderr);
         exit(1);
      } else {
         m=1;
         while(parm[m]!='\0') {
            switch (parm[m]) {
               case 'g':
                  ORTHO=ORTHO_GRAM;
                  break;
               case 'l':
                  ORTHO=ORTHO_LOWDIN;
                  break;
               case 'h':
                  DIAG=1;
                  break;
               default:
                  fprintf(stderr,"Switch '%c' unknown\n",parm[m]);
                  syntax(stderr);
                  exit(1);
                  break;
            }
            m++;
         }
      }
      k++;
   }
#ifdef DEBUG
   if(ORTHO==ORTHO_GRAM) fprintf(stderr,"Using Gram-Schmidt orthogonalization\n");
   else                  fprintf(stderr,"Using Lowdin orthogonalization\n");
#endif
#ifdef TRACE
   fprintf(stderr,"+++ Exiting parse\n");
#endif
   if(0) exit(1);
}
