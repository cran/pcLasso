// Automatically generated by SUtools, editing not advised.
#ifndef R_PCLASSO_H
#define R_PCLASSO_H
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("pclasso", String)
#else
#define _(String) (String)
#endif

#define FDEF(name)  {#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
void F77_SUB(pclasso)(
int *no,
int *ni,
double *x,
double *y,
double *w,
double *theta,
int *ng,
int *mg,
double *aa,
int *ne,
int *nx,
int *nlam,
double *ulam,
double *thr,
int *maxit,
int *verbose,
double *ao,
int *ia,
int *kin,
int *nlp,
int *jerr
);
 
static R_NativePrimitiveArgType pclasso_t[] = {
INTSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
INTSXP,
INTSXP,
INTSXP
};

void F77_SUB(logpclasso)(
int *no,
int *ni,
double *x,
double *y,
double *w,
double *theta,
int *ng,
int *mg,
double *aa,
int *ne,
int *nx,
int *nlam,
double *ulam,
double *thr,
int *maxit,
int *verbose,
double *a0,
double *ao,
int *ia,
int *kin,
int *nlp,
int *jerr
);
 
static R_NativePrimitiveArgType logpclasso_t[] = {
INTSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
INTSXP,
INTSXP
};

static R_FortranMethodDef fMethods[] = {
FDEF(pclasso) ,
FDEF(logpclasso) ,
{NULL, NULL, 0}
};

void R_init_pclasso(DllInfo *dll){
  R_registerRoutines(dll, NULL, NULL, fMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
#endif
