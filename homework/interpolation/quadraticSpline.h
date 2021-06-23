#ifndef HAVE_QUADSPLINE_H
#define HAVE_QUADSPLINE_H

#include <stdlib.h>
#include <assert.h>

typedef struct
{
    int numOfPts; /* n   */
    double *pts; /* x's */
    double *funcVals; /* y's */
    double *firstCoeff; /* b_i */
    double *secondCoeff; /* c_i */ } quadSpline;

quadSpline *initialize_quadratic_spline(int numberOfPoints, double *points, double *functionValues);

double evaluate_quadratic_spline(quadSpline *spline, double evaluationPoints);

double evaluate_quadratic_spline_derivative(quadSpline *spline, double evaluationPoints);

double evaluate_quadratic_spline_integral(quadSpline *spline, double evaluationPoints);

void free_quadratic_spline(quadSpline *spline);

#endif