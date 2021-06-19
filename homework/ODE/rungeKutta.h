#ifndef HAVE_RUNGEKUTTA_H
#define HAVE_RUNGEKUTTA_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void
rk_step12(void (*function)(double, gsl_vector *, gsl_vector *), double variable, gsl_vector *functionValue, double step,
          gsl_vector *functionStep, gsl_vector *error);

void rk_driver(void (*function)(double, gsl_vector *, gsl_vector *), double leftEndpoint, gsl_vector *functionValueLeft,
               double rightEndPoint, gsl_vector *functionValueRight, double step, double absoluteAccuracy,
               double relativeAccuracy, FILE *pathToFile);

#endif //HAVE_RUNGEKUTTA_H
