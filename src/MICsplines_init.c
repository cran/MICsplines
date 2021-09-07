#include <stdlib.h>  
#include <R_ext/Rdynload.h>

/* .C calls */
extern void MIC_splines_basis_C(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);


static const R_CMethodDef CEntries[] = {
    {"MIC_splines_basis_C", (DL_FUNC) &MIC_splines_basis_C, 15},
    {NULL, NULL, 0}
};

void R_init_MICsplines(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
