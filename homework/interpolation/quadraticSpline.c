#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "utilities.h"
#include "quadraticSpline.h"

quadSpline *initialize_quadratic_spline(int numberOfPoints, double *points, double *functionValueOfPoints)
{
    quadSpline *spline = (quadSpline *) malloc(sizeof(quadSpline));

    int numberOfEquations = numberOfPoints - 1;

    spline->numberOfPoints = numberOfPoints;
    spline->points = (double *) (numberOfEquations * sizeof(double));
    spline->functionValueOfPoints = (double *) (numberOfEquations * sizeof(double));
    spline->firstCoefficient = (double *) malloc(numberOfEquations * sizeof(double));
    spline->secondCoefficient = (double *) malloc(numberOfEquations * sizeof(double));

    for (int i = 0; i < numberOfPoints; i++)
    {
        spline->points[i] = points[i];
        spline->functionValueOfPoints[i] = functionValueOfPoints[i];
    }

    double pointsDifference[numberOfPoints - 1];
    double slope[numberOfPoints - 1];

    for (int i = 0; i < numberOfPoints - 1; i++)
    {
        pointsDifference[i] = points[i + 1] - points[i];
        slope[i] = (functionValueOfPoints[i + 1] - functionValueOfPoints[i]) / pointsDifference[i];
    }

    //Forward recursion to compute second coefficient
    spline->secondCoefficient[0] = 0;
    for (int i = 0; i < numberOfPoints - 2; i++)
    {
        spline->secondCoefficient[i + 1] =
                (slope[i + 1] - slope[i] - spline->secondCoefficient[i] * pointsDifference[i]) /
                pointsDifference[i + 1];
    }

    //Backward recursion to compute second coefficient
    spline->secondCoefficient[numberOfPoints - 2] /= 2;
    for (int i = numberOfPoints - 3; i >= 0; i--)
    {
        spline->secondCoefficient[i] =
                (slope[i + 1] - slope[i] - spline->secondCoefficient[i + 1] * pointsDifference[i + 1]) /
                pointsDifference[i];
    }

    //First coefficient
    for (int i = 0; i < numberOfPoints - 1; i++)
    {
        spline->firstCoefficient[i] = slope[i] - spline->secondCoefficient[i] * pointsDifference[i];
    }

    return spline;
}

double evaluate_quadratic_spline(quadSpline *spline, double pointsToEvaluateInterpolant)
{
    int whichInterval = binary_search(spline->numberOfPoints, spline->points, pointsToEvaluateInterpolant);
    double pointsDifference = pointsToEvaluateInterpolant - spline->points[whichInterval];
    double interpolantValue = spline->functionValueOfPoints[whichInterval] + pointsDifference *
                                                                             (spline->firstCoefficient[whichInterval] +
                                                                              pointsDifference *
                                                                              spline->secondCoefficient[whichInterval]);
    return interpolantValue;
}

double evaluate_quadratic_spline_derivative(quadSpline *spline, double pointsToEvaluateInterpolant)
{
    int whichInterval = binary_search(spline->numberOfPoints, spline->points, pointsToEvaluateInterpolant);
    double pointsDifference = pointsToEvaluateInterpolant - spline->points[whichInterval];
    double interpolantDerivativeValue =
            spline->firstCoefficient[whichInterval] + 2 * pointsDifference * spline->secondCoefficient[whichInterval];
    return interpolantDerivativeValue;
}

double evaluate_quadratic_spline_integral(quadSpline *spline, double pointsToEvaluateInterpolant)
{
    int whichInterval = binary_search(spline->numberOfPoints, spline->points, pointsToEvaluateInterpolant);
    double integralValue = 0;
    double pointsDifference;

    for (int i = 0; i <= whichInterval; i++)
    {
        if (i < whichInterval)
        {
            pointsDifference = spline->points[i + 1] - spline->points[i];
        }
        else
        {
            pointsDifference = pointsToEvaluateInterpolant - spline->points[i];
        }
        integralValue += spline->functionValueOfPoints[i] * pointsDifference +
                         spline->firstCoefficient[i] * pointsDifference * pointsDifference / 2 +
                         spline->secondCoefficient[i] * pointsDifference * pointsDifference * pointsDifference / 3;
    }
    return integralValue;
}

