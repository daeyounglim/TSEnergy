#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP boxcoxar_amh(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP calc_modelfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP singleMonitor(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// static const R_CallMethodDef CallEntries[] = {
//     {"boxcoxar_amh",  (DL_FUNC) &boxcoxar_amh,  16},
//     {"calc_modelfit", (DL_FUNC) &calc_modelfit, 12},
//     {"singleMonitor", (DL_FUNC) &singleMonitor, 10},
//     {NULL, NULL, 0}
// };

void R_init_TSEnergy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
