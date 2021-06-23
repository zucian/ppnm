#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "utilities.h"
#include "quadraticSpline.h"

quadSpline *initialize_quadratic_spline(int numberOfPoints, double *points, double *functionValues)
{
    quadSpline *spline = (quadSpline *) malloc(sizeof(quadSpline));

    int numberOfEquations = numberOfPoints - 1;

    spline->firstCoeff = (double *) malloc(numberOfEquations * sizeof(double));
    spline->secondCoeff = (double *) malloc(numberOfEquations * sizeof(double));
    spline->pts = (double *) malloc(numberOfPoints * sizeof(double));
    spline->funcVals = (double *) malloc(numberOfPoints * sizeof(double));
    spline->numOfPts = numberOfPoints;

    for (int i = 0; i < numberOfPoints; i++)
    {
        spline->pts[i] = points[i];
        spline->funcVals[i] = functionValues[i];
    }

    double pointsDifference[numberOfPoints - 1];
    double slope[numberOfPoints - 1];

    for (int i = 0; i < numberOfPoints - 1; i++)
    {
        pointsDifference[i] = points[i + 1] - points[i];
        slope[i] = (functionValues[i + 1] - functionValues[i]) / pointsDifference[i];
    }

    spline->secondCoeff[0] = 0;
    for (int i = 0; i < numberOfPoints - 2; i++)
    {
        spline->secondCoeff[i + 1] =
                (slope[i + 1] - slope[i] - (spline->secondCoeff[i]) * pointsDifference[i]) / pointsDifference[i + 1];
    }

    spline->secondCoeff[numberOfPoints - 2] /= 2;
    for (int i = numberOfPoints - 3; i >= 0; i--)
    {
        spline->secondCoeff[i] =
                (slope[i + 1] - slope[i] - (spline->secondCoeff[i + 1]) * pointsDifference[i + 1]) /
                pointsDifference[i];
    }

    for (int i = 0; i < numberOfPoints - 1; i++)
    {
        spline->firstCoeff[i] = slope[i] - (spline->secondCoeff[i]) * pointsDifference[i];
    }

    return spline;
}

double evaluate_quadratic_spline(quadSpline *spline, double evaluationPoints)
{
    assert((evaluationPoints >= (spline->pts[0])) && (evaluationPoints <= (spline->pts[spline->numOfPts - 1])));

    int whichInterval = binarySearch(spline->numOfPts, spline->pts, evaluationPoints);
    double pointsDifference = evaluationPoints - (spline->pts[whichInterval]);

    double interpolantValue = (spline->funcVals[whichInterval]) +
                              pointsDifference * ((spline->firstCoeff[whichInterval]) +
                                                  pointsDifference * (spline->secondCoeff[whichInterval]));
    return interpolantValue;
}


double evaluate_quadratic_spline_derivative(quadSpline *spline, double evaluationPoints)
{
    assert((evaluationPoints >= (spline->pts[0])) && (evaluationPoints <= (spline->pts[spline->numOfPts - 1])));

    //  Use a binary search to determine which subinterval evalPt is in
    int whichInterval = binarySearch(spline->numOfPts, spline->pts, evaluationPoints);
    double pointsDifference = evaluationPoints - (spline->pts[whichInterval]);

    double interpolantValue =
            (spline->firstCoeff[whichInterval]) + 2 * pointsDifference * (spline->secondCoeff[whichInterval]);
    return interpolantValue;
}

double evaluate_quadratic_spline_integral(quadSpline *spline, double evaluationPoints)
{
    int whichInterval = binarySearch((spline->numOfPts), (spline->pts), evaluationPoints);

    double integral = 0;
    double pointsDifference;
    for (int i = 0; i <= whichInterval; i++)
    {

        // Note: pow(x, 2) is slower than x*x in principle, but they are definitively comparable below 1e6 elements
        if (i < whichInterval)
        {
            pointsDifference = ((spline->pts[i + 1]) - (spline->pts[i]));
        }
        else
        {
            pointsDifference = (evaluationPoints - (spline->pts[i]));
        }
        integral += (spline->funcVals[i]) * pointsDifference +
                    (spline->firstCoeff[i]) * pointsDifference * pointsDifference / 2 +
                    (spline->secondCoeff[i]) * pointsDifference * pointsDifference * pointsDifference / 3;
    }

    return integral;
}


void free_quadratic_spline(quadSpline *spline)
{
    free(spline->pts);
    free(spline->funcVals);
    free(spline->firstCoeff);
    free(spline->secondCoeff);
    free(spline);
}