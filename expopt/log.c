/****************************************************************************/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include "expopt.h"

static const char *shell_sym[] = {"s","p","d","f","g","h","i","k","l","m"};
static const char *gen_name[]  = {"?","HotTempered","Harmonic","Polynomial",
                                   "Cosine","Unit","Exponential","Input",
                                   "SplitValence","HotPoly"};
static const char *meth_name[] = {"?","Auto","Amoeba","Powell","QN","QN2"};

void OptimizationLog(int opt_idx) {
   char fname[32];
   int l = SYS.active;
   snprintf(fname, sizeof(fname), "opt_%03i.log", opt_idx);
   SYS.logfile = fopen(fname, "w");
   if(SYS.logfile == NULL) {
      fprintf(stderr, "Could not open log file %s\n", fname);
      SYS.logfile = stdout;
   }
   fprintf(SYS.logfile,
      "Opt[%i] shell %i(%s) Z=%i n=%i dirs=%i  gen=%s  method=%s%s\n",
      opt_idx, l, shell_sym[l], SYS.basis[l].Z,
      SYS.basis[l].n, SYS.n_dir,
      gen_name[SYS.generator], meth_name[SYS.method],
      SYS.no_linesearch ? " [NoLineSearch]" : "");
   if(SYS.no_linesearch)
      fprintf(SYS.logfile, "%4s  %6s  %8s  %14s  %12s  %10s  %8s  %8s\n",
         "iter","fevals","wall(s)","energy","dE","gnorm","lambda","trust");
   else
      fprintf(SYS.logfile, "%4s  %6s  %6s  %8s  %14s  %12s  %10s  %8s\n",
         "iter","fevals","lsstp","wall(s)","energy","dE","gnorm","lambda");
}
