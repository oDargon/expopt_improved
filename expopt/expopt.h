#ifndef __EXPOPT_H__
#define __EXPOPT_H__
#define VERSION "2.2a - VERSION UNDER DEVELOPMENT!!!"
/*--------------------------------------------------------------------------*/
/* Limits                                                                   */
/*--------------------------------------------------------------------------*/
#define MAX_DIR   32
#define MAX_LQN    9
/*--------------------------------------------------------------------------*/
/* Switches                                                                 */
/*--------------------------------------------------------------------------*/
#define GEN_HOT       1
#define GEN_HARM      2
#define GEN_POLY      3
#define GEN_COS       4
#define GEN_UNIT      5
#define GEN_EXP       6
#define GEN_INPUT     7
#define GEN_SPLIT     8
#define GEN_HOTPO     9
#define OPT_AUTO      1
#define OPT_AMOEBA    2
#define OPT_POWELL    3
#define OPT_QN        4
#define OPT_QN2       5
/*--------------------------------------------------------------------------*/
/* Typedefs                                                                 */
/*--------------------------------------------------------------------------*/
typedef struct {
   int      Z,l,n;
   double  *z;
} shell;
typedef struct {
   double  charge;
   double  ethr;
   double  gthr;
   double  energy;
   double  pen_a;
   double  pen_b;
   double  pen_c;
   double  B0;
   double  A0;
   double  alpha;
   int     fnc_eval;
   int     verbose;
   int     generator;
   int     method;
   int     active;
   int     n_dir;
   int    *index;
   int     printthis;
   int     printraw;
   shell   basis[MAX_LQN+1];
   double  directions[MAX_DIR*MAX_DIR];
} sysdef;
/*--------------------------------------------------------------------------*/
/* Global variables.                                                        */
/*--------------------------------------------------------------------------*/
extern sysdef    SYS;
/*--------------------------------------------------------------------------*/
/* external functions.                                                      */
/*--------------------------------------------------------------------------*/
extern void    input_control(FILE *, FILE *);
extern double  energy(double *);
extern void    optimize(int, int *, int, FILE *);
extern void    amoeba(double *, int, double, FILE *);
extern void    quasi_newton(double *, int, double, double, FILE *);
/* extern void    qn2(double *, int, double, double, FILE *); */
extern void    powell(double *, int, double, FILE *);
extern double  make_scale(int, double *);
extern void    gauss(double *, double *, double *, int);
extern double  line_search(double *, double *, int n, double *, double, FILE *);
extern int     find_min(double *, int);
extern void    guessian(double *, double *, int, FILE *);
extern void    extend(int, int, FILE *);
extern void    check_hessian(double *, int);
extern void    jacobi(double *, double *, int, int);
extern double  check_lin_dep(double M[MAX_DIR][MAX_DIR], int , FILE *);
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
#endif
