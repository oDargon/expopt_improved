/**************************************************************************/
/*                                                                        */
/*                                                                        */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/*                                                                        */
/*                                                                        */
/**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "onbas.h"
int ORTHO;
int DIAG;
/*========================================================================*/
/*                                                                        */
/*========================================================================*/
int main(int argc, char *argv[]) {
   basis   X;
#ifdef TRACE
   fprintf(stderr,"+++ Entering main\n");
#endif
   parse(argc,argv);
   rdbas(stdin,&X);
#ifdef DEBUG
   {
      int lqn;
      fprintf(stderr,"Z=%.0f, maxlqn=%i\n",X.Z,X.maxlqn);
      for(lqn=0; lqn<=X.maxlqn; lqn++) {
         fprintf(stderr,"nprim[%i]=%i, ncont[%i]=%i\n",lqn,X.shell[lqn].nprim,lqn,X.shell[lqn].ncont);
      }
   }
#endif
   if(ORTHO==ORTHO_GRAM) {
      gram(&X);
   } else if(ORTHO==ORTHO_LOWDIN) {
      lowdin(&X);
   } else {
      fprintf(stderr,"Internal inconsistency: ORTHO=%i\n",ORTHO);
      exit(1);
   }
   if(DIAG) diag(&X);
   prbas(stdout,X);
#ifdef TRACE
   fprintf(stderr,"+++ Exiting main\n");
#endif
   return 0;
}
