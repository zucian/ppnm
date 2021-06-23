#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "utilities.h"
#include "cubicSpline.h"


cubicSpline *initialize_cubic_spline(int numberOfPoints, double *points, double *functionValues)
{
    int numberOfEquations = numberOfPoints - 1;
    cubicSpline *spline = (cubicSpline *) malloc(sizeof(cubicSpline));
    spline->pts = (double *) malloc(numberOfPoints * sizeof(double));
    spline->funcVals = (double *) malloc(numberOfPoints * sizeof(double));
    spline->firstCoeff = (double *) malloc(numberOfPoints * sizeof(double));
    spline->secondCoeff = (double *) malloc(numberOfEquations * sizeof(double));
    spline->thirdCoeff = (double *) malloc(numberOfEquations * sizeof(double));
    spline->numOfPts = numberOfPoints;

    // Fill field values from input data
    for (int i = 0; i < numberOfPoints; i++)
    {
        spline->pts[i] = points[i];
        spline->funcVals[i] = functionValues[i];
    }

    // Set up p_i h_i variables
    double pointsDifference[numberOfEquations];
    double slope[numberOfEquations];
    for (int i = 0; i < numberOfEquations; i++)
    {
        pointsDifference[i] = points[i + 1] - points[i];
        assert(pointsDifference[i] > 0);

        slope[i] = (functionValues[i + 1] - functionValues[i]) / pointsDifference[i];
    }

    //set up system of equation (upper triangular system)
    double matrixDiagonal[numberOfPoints];
    double matrixSuperDiagonal[numberOfPoints - 1];
    double resultVector[numberOfPoints];

    //Recusrion start
    matrixDiagonal[0] = 2;
    matrixSuperDiagonal[0] = 1;
    resultVector[0] = 3 * slope[0];

    //array filled according to recursion
    for (int i = 0; i < numberOfEquations - 1; i++)
    {
        matrixDiagonal[i + 1] = 2 * pointsDifference[i] / pointsDifference[i + 1] + 2;
        matrixSuperDiagonal[i + 1] = pointsDifference[i] / pointsDifference[i + 1];
        resultVector[i + 1] = 3 * (slope[i] + slope[i + 1] * pointsDifference[i] / pointsDifference[i + 1]);
    }
    matrixDiagonal[numberOfPoints - 1] = 2;
    resultVector[numberOfPoints - 1] = 3 * slope[numberOfPoints - 2];

    //Gaus elimination eq 29
    for (int i = 1; i < numberOfPoints; i++)
    {
        matrixDiagonal[i] -= matrixSuperDiagonal[i - 1] / matrixDiagonal[i - 1];
        resultVector[i] -= resultVector[i - 1] / matrixDiagonal[i - 1];
    }

    //Backward substitution for b_i
    spline->firstCoeff[numberOfEquations] = resultVector[numberOfPoints - 1] / matrixDiagonal[numberOfPoints - 1];
    for (int i = numberOfEquations - 1; i >= 0; i--)
    {
        spline->firstCoeff[i] =
                (resultVector[i] - matrixSuperDiagonal[i] * (spline->firstCoeff[i + 1])) / matrixDiagonal[i];
    }

    //Compute remaining using eq (20)
    for (int i = 0; i < numberOfEquations; i++)
    {
        spline->secondCoeff[i] =
                (-2 * (spline->firstCoeff[i]) - (spline->firstCoeff[i + 1]) + 3 * slope[i]) / pointsDifference[i];
        spline->thirdCoeff[i] =
                ((spline->firstCoeff[i]) + (spline->firstCoeff[i + 1]) - 2 * slope[i]) / pointsDifference[i] /
                pointsDifference[i];
    }

    return spline;
}

double evaluate_cubic_spline(cubicSpline *spline, double evaluationPoint)
{
    assert((evaluationPoint >= (spline->pts[0])) && (evaluationPoint <= (spline->pts[spline->numOfPts - 1])));

    //Binary search to find subinterval for evaluation point
    int whichInterval = binarySearch(spline->numOfPts, spline->pts, evaluationPoint);

    double pointsDifference = evaluationPoint - (spline->pts[whichInterval]);
    double thirdDifference = pointsDifference * (spline->thirdCoeff[whichInterval]);
    double secondDifference = pointsDifference * ((spline->secondCoeff[whichInterval]) + thirdDifference);
    double firstDifference = pointsDifference * ((spline->firstCoeff[whichInterval]) + secondDifference);

    double interpolantValue = (spline->funcVals[whichInterval]) + firstDifference;
    return interpolantValue;
}

double evaluate_cubic_spline_derivative(cubicSpline *spline, double evaluationPoint)
{
    assert((evaluationPoint >= (spline->pts[0])) && (evaluationPoint <= (spline->pts[spline->numOfPts - 1])));

    //  Use a binary search to determine which subinterval evalPt is in
    int whichInterval = binarySearch(spline->numOfPts, spline->pts, evaluationPoint);
    double pointsDifference = evaluationPoint - (spline->pts[whichInterval]);

    double interpolantValue =
            (spline->firstCoeff[whichInterval]) + 2 * pointsDifference * (spline->secondCoeff[whichInterval]) +
            3 * pointsDifference * pointsDifference * (spline->thirdCoeff[whichInterval]);
    return interpolantValue;
}

double evaluate_cubic_spline_integral(cubicSpline *spline, double evaluationPoint)
{

    // Binary search to find subinterval for evaluation point
    int whichInterval = binarySearch((spline->numOfPts), (spline->pts), evaluationPoint);

    double integral = 0;
    double pointsDifference;
    for (int i = 0; i <= whichInterval; i++)
    {

        if (i < whichInterval)
        {
            pointsDifference = ((spline->pts[i + 1]) - (spline->pts[i]));
        }
        else
        {
            pointsDifference = (evaluationPoint - (spline->pts[i]));
        }
        integral += (spline->funcVals[i]) * pointsDifference +
                    (spline->firstCoeff[i]) * pointsDifference * pointsDifference / 2 +
                    (spline->secondCoeff[i]) * pointsDifference * pointsDifference * pointsDifference / 3 +
                    (spline->thirdCoeff[i]) * pointsDifference * pointsDifference * pointsDifference *
                    pointsDifference / 4;
    }

    return integral;
}


void free_cubic_spline(cubicSpline *spline)
{
    free(spline->pts);
    free(spline->funcVals);
    free(spline->firstCoeff);
    free(spline->secondCoeff);
    free(spline->thirdCoeff);
    free(spline);
}