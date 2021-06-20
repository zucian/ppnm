#ifndef HAVE_ROOTFINDINGMETHODS_H
#define HAVE_ROOTFINDINGMETHODS_H

#include <gsl/gsl_vector.h>

void
newton_method(void function(gsl_vector *point, gsl_vector *functionValues), gsl_vector *startingPoint, double tolerance);

#endif //HAVE_ROOTFINDINGMETHODS_H
