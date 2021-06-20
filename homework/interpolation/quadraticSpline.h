#ifndef HAVE_QUADRATICSPLINE_H
#define HAVE_QUADRATICSPLINE_H

#include <stdlib.h>
#include <assert.h>

typedef struct
{
    int numberOfPoints;
    double *points;
    double *functionValueOfPoints;
    double *firstCoefficient;
    double *secondCoefficient;
} quadSpline;

quadSpline *initialize_quadratic_spline(int numberOfPoints, double *points, double *functionValueOfPoints);

double evaluate_quadratic_spline(quadSpline *spline, double pointsToEvaluateInterpolant);

double evaluate_quadratic_spline_derivative(quadSpline *spline, double pointsToEvaluateInterpolant);

double evaluate_quadratic_spline_integral(quadSpline *spline, double pointsToEvaluateInterpolant);

void free_quadratic_spline(quadSpline *spline);

#endif //HAVE_QUADRATICSPLINE_H