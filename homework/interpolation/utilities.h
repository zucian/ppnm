#ifndef HAVE_UTILITIES_H
#define HAVE_UTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>

int binary_search(int numOfPts, double *pts, double evalPt);

double integrate_GSL_function(double lowerLimit, double upperLimit, const gsl_function* gslFunction, double toleranceAbsolute, double toleranceRelative, size_t iterationLimit);

void input_to_array(double *XData, double *YData, char *inputFilename);

#endif