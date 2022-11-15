/**************************************************************************/
/*                                                                        */
/*                              J A C O R D                               */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* This subroutine sorts eigen-vectors with respect to eigen-values.      */
/*                                                                        */
/**************************************************************************/
void jacord(double f[], double u[], int n, int nr) {
   double t;
   int i,j,k;
   for(i=0; i<n; i++) {
      j=i;
      for(k=i; k<n; k++) if(f[k]<f[j]) j=k;
      t    = f[i];
      f[i] = f[j];
      f[j] = t;
      for(k=0; k<nr; k++) {
         t         =   u[k+i*nr];
         u[k+i*nr] = - u[k+j*nr];
         u[k+j*nr] =   t;
      }
   }
}
