/**************************************************************************/
/*                                                                        */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/*                                                                        */
/**************************************************************************/
#ifndef __ONBAS_H__
#define __ONBAS_H__
/*========================================================================*/
/* Defines                                                                */
/*========================================================================*/
#define MXPRIM  32
#define MXCONT  32
#define MXLQN    6
#define ORTHO_GRAM    1
#define ORTHO_LOWDIN  2
#define PI 3.14159265359
/*========================================================================*/
/* Type definitions                                                       */
/*========================================================================*/
typedef struct {
   double z[MXPRIM];
   double c[MXPRIM][MXCONT];
   int    nprim;
   int    ncont;
} shelltype;
typedef struct {
   shelltype shell[MXLQN+1];
   double    Z;
   int       En;
   int       maxlqn;
   double    Enc[MXPRIM];
} basis;
/*========================================================================*/
/* Globals                                                                */
/*========================================================================*/
extern int ORTHO;
extern int DIAG;
/*========================================================================*/
/* Externals                                                              */
/*========================================================================*/
extern void parse(int argc, char *argv[]);
extern void syntax(FILE *f);
extern void rdbas(FILE *f, basis *Y);
extern void prbas(FILE *f, basis X);
extern void gram(basis *Y);
extern void lowdin(basis *Y);
extern void diag(basis *Y);
extern void jacobi(double f[], double u[], int n, int nr);
extern void jacord(double f[], double u[], int n, int nr);
/*========================================================================*/
/*                                                                        */
/*========================================================================*/
#endif
