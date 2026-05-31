/****************************************************************************/
/*                                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "expopt.h"
/*--------------------------------------------------------------------------*/
typedef struct {
   int    active;   /* 0=free  1=prepared  2=done */
   pid_t  pid;
   char   dir[4096];
   int    n_s;
   double *s;
} JobRecord;
static JobRecord *job_table  = NULL;
static int        jobs_init   = 0;
static char       run_dir[4096-64];
/*--------------------------------------------------------------------------*/
static void init_jobs(void) {
   int k;
   int table_size=SYS.max_jobs > 2*MAX_DIR+2 ? SYS.max_jobs : 2*MAX_DIR+2;
   job_table=(JobRecord *)malloc(table_size*sizeof(JobRecord));
   if(job_table==NULL) {
      fprintf(stderr,"Could not allocate job table\n");
      exit(1);
   }
   for(k=0; k<table_size; k++) job_table[k].active=0;
   if(getcwd(run_dir,sizeof(run_dir))==NULL) {
      fprintf(stderr,"Could not get working directory\n");
      exit(1);
   }
   char jobs_dir[4096];
   snprintf(jobs_dir,sizeof(jobs_dir),"%s/expopt_jobs",run_dir);
   mkdir(jobs_dir,0755);
   jobs_init=1;
}
/*--------------------------------------------------------------------------*/
static void write_basis(FILE *f, double *s) {
   int    k,l,m,n,max_lqn;
   double z;
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
               if(SYS.active==l) z=z*make_scale(k,s);
               fprintf(f,"%.8f\n",z);
            }
            if(frozen) {
               for(m=0; m<n; m++) {
                  for(k=0; k<nc; k++) fprintf(f,"%.8f ",SYS.basis[l].c[k*n+m]);
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
   } else {
      l=SYS.printthis;
      n=SYS.basis[l].n;
      if(SYS.printraw) {
         for(k=0; k<n; k++) {
            z=SYS.basis[l].z[k];
            if(SYS.active==l) z=z*make_scale(k,s);
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
               if(SYS.active==l) z=z*make_scale(k,s);
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
}
/*==========================================================================*/
int energy_prepare(double *s) {
   FILE  *f;
   char   basis_path[4128];
   char   cp_cmd[8224];
   int    slot,k;
   if(!jobs_init) init_jobs();
   int table_size=SYS.max_jobs > 2*MAX_DIR+2 ? SYS.max_jobs : 2*MAX_DIR+2;
   for(slot=0; slot<table_size; slot++)
      if(!job_table[slot].active) break;
   if(slot==table_size) {
      fprintf(stderr,"All job slots occupied — call energy_run first\n");
      exit(1);
   }
   snprintf(job_table[slot].dir,sizeof(job_table[slot].dir),"%s/expopt_jobs/expopt_job_%04i",run_dir,slot);
   mkdir(job_table[slot].dir,0755);
   snprintf(cp_cmd,sizeof(cp_cmd),"cp %s/energy.ksh %s/",run_dir,job_table[slot].dir);
   if(system(cp_cmd)!=0) {
      fprintf(stderr,"Could not copy energy.ksh to job directory\n");
      exit(1);
   }
   snprintf(basis_path,sizeof(basis_path),"%s/BASIS.DAT",job_table[slot].dir);
   if((f=fopen(basis_path,"w"))==NULL) {
      fprintf(stderr,"Could not open %s\n",basis_path);
      exit(1);
   }
   write_basis(f,s);
   fclose(f);
   if(s!=NULL) {
      job_table[slot].n_s=SYS.n_dir;
      job_table[slot].s=(double *)malloc(SYS.n_dir*sizeof(double));
      for(k=0; k<SYS.n_dir; k++) job_table[slot].s[k]=s[k];
   } else {
      job_table[slot].n_s=0;
      job_table[slot].s=NULL;
   }
   job_table[slot].active=1;
   return slot;
}
/*==========================================================================*/
void energy_run(int *handles, int n) {
   char   cmd[4128];
   int    status,k,h;
   int    queue_pos,running;
   pid_t  pid;
   fflush(stdout);
   if(SYS.max_jobs==1) {
      for(k=0; k<n; k++) {
         h=handles[k];
         snprintf(cmd,sizeof(cmd),"cd '%s' && ./energy.ksh",job_table[h].dir);
         if(system(cmd)!=0) {
            fprintf(stderr,"Error executing energy.ksh\n");
            exit(1);
         }
         job_table[h].active=2;
      }
   } else {
      queue_pos=0;
      running  =0;
      while(queue_pos<n || running>0) {
         while(running<SYS.max_jobs && queue_pos<n) {
            h=handles[queue_pos++];
            pid=fork();
            if(pid<0) {
               fprintf(stderr,"fork() failed\n");
               exit(1);
            }
            if(pid==0) {
               chdir(job_table[h].dir);
               execl("/bin/sh","sh","./energy.ksh",NULL);
               exit(1);
            }
            job_table[h].pid=pid;
            running++;
         }
         pid=waitpid(-1,&status,0);
         if(!WIFEXITED(status)||WEXITSTATUS(status)!=0) {
            fprintf(stderr,"Error executing energy.ksh\n");
            exit(1);
         }
         int ts=SYS.max_jobs > 2*MAX_DIR+2 ? SYS.max_jobs : 2*MAX_DIR+2;
         for(k=0; k<ts; k++) {
            if(job_table[k].active==1 && job_table[k].pid==pid) {
               job_table[k].active=2;
               break;
            }
         }
         running--;
      }
   }
}
/*==========================================================================*/
double energy_collect(int handle) {
   FILE   *f;
   char    line[128];
   char    energy_path[4128];
   int     k,l,n;
   double  e,z,p;
   snprintf(energy_path,sizeof(energy_path),"%s/ENERGY.DAT",job_table[handle].dir);
   if((f=fopen(energy_path,"r"))==NULL) {
      fprintf(stderr,"Could not open file '%s'\n",energy_path);
      exit(1);
   }
   if(fgets(line,sizeof(line),f)==NULL) {
      fprintf(stderr,"Premature end of ENERGY.DAT\n");
      exit(1);
   }
   e=atof(line);
   fclose(f);
   if(SYS.active>-1) {
      p=0.0;
      l=SYS.active;
      n=SYS.basis[l].n;
      for(k=0; k<n; k++) {
         z=SYS.basis[l].z[k];
         z=z*make_scale(k,job_table[handle].s);
         z=log(z);
         p=p+SYS.pen_a*exp(SYS.pen_b*(z-SYS.pen_c));
      }
      e=e+p;
   }
   if((SYS.active>-1)&&(e<SYS.energy)) {
      l=SYS.active;
      n=SYS.basis[l].n;
      printf("Energy %.8f\n",e);
      printf(" %2i ",SYS.basis[l].Z);
      printf("%2i ",SYS.basis[l].l);
      printf("%2i  ",SYS.basis[l].n);
      for(k=0; k<n; k++) {
         if((k!=0)&&(k%5==0)) printf("\n           ");
         z=SYS.basis[l].z[k];
         z=z*make_scale(k,job_table[handle].s);
         printf("%.8f ",z);
      }
      printf("\n");
      SYS.energy=e;
   }
   SYS.fnc_eval++;
   if(job_table[handle].s!=NULL) free(job_table[handle].s);
   job_table[handle].active=0;
   return e;
}
/*==========================================================================*/
double energy(double *s) {
   int handle=energy_prepare(s);
   energy_run(&handle,1);
   return energy_collect(handle);
}
