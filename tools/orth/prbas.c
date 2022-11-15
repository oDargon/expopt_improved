/****************************************************************************/
/*                                                                          */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include "onbas.h"
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void prbas(FILE *f, basis X) {
   int lqn;
   int ip,ic;
#ifdef TRACE
   fprintf(stderr,"+++ Entering prbas\n");
#endif
   fprintf(f," %.1f %i\n",X.Z,X.maxlqn);
   for(lqn=0; lqn<=X.maxlqn; lqn++) {
      fprintf(f," %i %i\n",X.shell[lqn].nprim,X.shell[lqn].ncont);
      for(ip=0; ip<X.shell[lqn].nprim; ip++) {
         if( ((ip%5)==0) && (ip!=0) ) fprintf(f,"\n");
         fprintf(f," %.8f",X.shell[lqn].z[ip]);
      }
      fprintf(f,"\n");
      for(ip=0; ip<X.shell[lqn].nprim; ip++) {
         for(ic=0; ic<X.shell[lqn].ncont; ic++) fprintf(f," %15.10f",X.shell[lqn].c[ip][ic]);
         fprintf(f,"\n");
      }
   }
#ifdef TRACE
   fprintf(stderr,"+++ Exiting prbas\n");
#endif
}
