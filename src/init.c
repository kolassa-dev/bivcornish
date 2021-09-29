#include <R_ext/RS.h>                                                                                                 
#include <stdlib.h> // for NULL                                                                                       
#include <R_ext/Rdynload.h>                                                                                           

/* FIXME:                                                                                                             
   Check these declarations against the C/Fortran source code.                                                        
*/                                                                                      

/* .Fortran calls */
extern void F77_NAME(bivcorn)(void *, void *, void *, void *, void *, void *, 
void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
void *, void *, void *);
extern void F77_NAME(bivtail)(void *, void *, void *, void *, void *, void *, 
void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(isecnorm)(void *, void *, void *, void *, void *);                                               

static const R_FortranMethodDef FortranEntries[] = {                                                                  
    {"bivcorn",  (DL_FUNC) &F77_NAME(bivcorn),  19},                                                                  
    {"bivtail",  (DL_FUNC) &F77_NAME(bivtail),  14},                                                                  
    {"isecnorm", (DL_FUNC) &F77_NAME(isecnorm),  5},
    {NULL, NULL, 0}
};

void R_init_bivcornish(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
