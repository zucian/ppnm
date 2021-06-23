#ifndef HAVE_LINEARSPLINE_H
#define HAVE_LINEARSPLINE_H

double linear_spline(int numberOfPoints, double *points, double *functionValues, double evaluationPoints);

double linear_spline_integrate(int numberOfPoints, double *points, double *functionValues, double evaluationPoint);

#endif