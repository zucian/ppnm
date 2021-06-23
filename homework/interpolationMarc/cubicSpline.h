#ifndef HAVE_CUBICSPLINE_H
#define HAVE_CUBICSPLINE_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct
{
    int numOfPts; /* n   */
    double *pts; /* x's */
    double *funcVals; /* y's */
    double *firstCoeff; /* b_i */
    double *secondCoeff; /* c_i */
    double *thirdCoeff; /* d_i */ } cubicSpline;

cubicSpline *initialize_cubic_spline(int numberOfPoints, double *points, double *functionValues);

double evaluate_cubic_spline(cubicSpline *spline, double evaluationPoint);

double evaluate_cubic_spline_derivative(cubicSpline *spline, double evaluationPoint);

double evaluate_cubic_spline_integral(cubicSpline *spline, double evaluationPoint);

void free_cubic_spline(cubicSpline *spline);

#endif