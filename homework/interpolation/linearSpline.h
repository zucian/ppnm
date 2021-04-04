#ifndef HAVE_LINEARSPLINE_H
#define HAVE_LINEARSPLINE_H

double linear_spline(int numOfPts, double *pts, double *funcVals, double evalPt);

double linear_spline_integration(int numOfPts, double *pts, double *funcVals, double evalPt);

#endif