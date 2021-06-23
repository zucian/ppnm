#ifndef HAVE_LINSPLINE_H
#define HAVE_LINSPLINE_H

double linear_spline(int numberOfPoints, double *points, double *functionValues, double evaluationPoints);

double linear_spline_integrate(int numberOfPoints, double *points, double *functionValues, double evaluationPoint);

#endif