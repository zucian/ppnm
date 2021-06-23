#ifndef HAVE_GSLINTFUNC_H
#define HAVE_GSLINTFUNC_H

#include <gsl/gsl_integration.h>

double integrate_function( double lowerLimit, double upperLimit, const gsl_function* myGSLFUNC, double tolAbs, double tolRel, size_t limit );

#endif