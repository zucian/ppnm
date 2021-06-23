#ifndef HAVE_CUBICSPLINE_H
#define HAVE_CUBICSPLINE_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct
{
    int numberOfPoints;
    double *points;
    double *functionValues;
    double *firstCoefficient;
    double *secondCoefficient;
    double *thirdCoefficient;
} cubicSpline;

cubicSpline *initialize_cubic_spline(int numberOfPoints, double *points, double *functionValues);

double evaluate_cubic_spline(cubicSpline *spline, double evaluationPoint);

double evaluate_cubic_spline_derivative(cubicSpline *spline, double evaluationPoint);

double evaluate_cubic_spline_integral(cubicSpline *spline, double evaluationPoints);

void free_cubic_spline(cubicSpline *spline);

#endif //HAVE_CUBICSPLINE_H
