#ifndef UNI_SUBSPLINE_H
#define UNI_SUBSPLINE_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct
{
    int numberOfPoints; /* n  */
    double *points; /* x's */
    double *functionValues; /* y's */
    double *derivativeValues; /* p's */
    double *firstCoefficient; /* b_i */
    double *secondCoefficient; /* c_i */
    double *thirdCoefficient; /* d_i */ } subSpline;

subSpline *initialize_sub_spline(int numberOfPoints, double *xData, double *yData, double *pData);

double evaluate_sub_spline(subSpline *spline, double evaluationPoint);

void free_sub_spline(subSpline *spline);

#endif //UNI_SUBSPLINE_H
