#ifndef HAVE_INPUTARRAY_H
#define HAVE_INPUTARRAY_H

#include <gsl/gsl_vector.h>

void input_to_array(int numberOfDataPoints, gsl_vector *xData, gsl_vector *yData, char *inputFilename);

#endif
