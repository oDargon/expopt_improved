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
void syntax(FILE *f) {
#ifdef TRACE
   fprintf(stderr,"+++ Entering syntax\n");
#endif
   fprintf(f,"\n");
   fprintf(f,"Syntax: onbas <options> \n");
   fprintf(f,"  options: -g  -- Use Gram-Schmidt orthogonalization\n");
   fprintf(f,"           -l  -- Use Lowdin orthogonalization\n");
   fprintf(f,"           -h  -- Diagonalize one electron hamiltonian\n");
   fprintf(f,"\n");
#ifdef TRACE
   fprintf(stderr,"+++ Exiting syntax\n");
#endif
}
