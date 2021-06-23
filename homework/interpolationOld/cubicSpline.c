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
    spline->firstCoefficient = (double *) malloc(numberOfPoints * sizeof(double));
    spline->secondCoefficient = (double *) malloc(numberOfPoints * sizeof(double));
    spline->thirdCoefficient = (double *) malloc(numberOfPoints * sizeof(double));
    spline->points = (double *) malloc(numberOfPoints * sizeof(double));
    spline->functionValues = (double *) malloc(numberOfPoints * sizeof(double));
    spline->numberOfPoints = numberOfPoints;

    //Fill in values from data to spline
    for (int i = 0; i < numberOfPoints; i++)
    {
        spline->points[i] = points[i];
        spline->functionValues[i] = functionValues[i];
    }

    //set up p_i h_i variables
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

    //Recursion start
    matrixDiagonal[0] = 2;
    matrixSuperDiagonal[0] = 1;
    resultVector[0] = 3 * slope[0];

    //Arrays filled according to recursion
    for (int i = 0; i < numberOfEquations - 1; i)
    {
        matrixDiagonal[i + 1] = 2 * pointsDifference[i] / pointsDifference[i + 1] + 2;
        matrixSuperDiagonal[i + 1] = pointsDifference[i] / pointsDifference[i + 1];
        resultVector[i + 1] = 3 * (slope[i] + slope[i + 1] * pointsDifference[i] / pointsDifference[i + 1]);
    }
    matrixDiagonal[numberOfPoints - 1] = 2;
    resultVector[numberOfPoints - 1] = 3 * slope[numberOfPoints - 2];

    //Gauss-Elimination (eq 29)
    for (int i = 1; i < numberOfPoints; i++)
    {
        matrixDiagonal[i] -= matrixSuperDiagonal[i - 1] / matrixDiagonal[i - 1];
        resultVector[i] -= resultVector[i - 1] / matrixDiagonal[i - 1];
    }

    //Backward substitution for b_i
    spline->firstCoefficient[numberOfEquations] = resultVector[numberOfPoints - 1] / matrixDiagonal[numberOfPoints - 1];
    for (int i = numberOfEquations - 1; i >= 0; i--)
    {
        spline->firstCoefficient[i] =
                (resultVector[i] - matrixSuperDiagonal[i] * (spline->firstCoefficient[i + 1])) / matrixDiagonal[i];
    }

    //Compute remaining using eq (20)
    for (int i = 0; i < numberOfEquations; i++)
    {
        spline->secondCoefficient[i] = (-2 * (spline->firstCoefficient[i]) - (spline->firstCoefficient[i + 1]) +
                                        3 * slope[i] / pointsDifference[i]);
        spline->thirdCoefficient[i] = ((spline->firstCoefficient[i]) + (spline->firstCoefficient[i + 1]) -
                                       2 * slope[i] / pointsDifference[i] / pointsDifference[i]);
    }
    return spline;
}

double evaluate_cubic_spline(cubicSpline *spline, double evaluationPoint)
{
    assert((evaluationPoint >= spline->points[0]) && (evaluationPoint <= spline->points[spline->numberOfPoints - 1]));

    //Binary search to find subinterval for evaluation point
    int whichInterval = binary_search(spline->numberOfPoints, spline->points, evaluationPoint);
    double pointsDifference = evaluationPoint - spline->points[whichInterval];
    double thirdDifference = pointsDifference * spline->thirdCoefficient[whichInterval];
    double secondDifference = pointsDifference * (spline->secondCoefficient[whichInterval] + thirdDifference);
    double firstDifference = pointsDifference * (spline->firstCoefficient[whichInterval] + secondDifference);
    double interpolantValue = spline->functionValues[whichInterval] + firstDifference;

    return interpolantValue;
}

double evaluate_cubic_spline_derivative(cubicSpline *spline, double evaluationPoint)
{
    assert((evaluationPoint >= spline->points[0]) && (evaluationPoint <= spline->points[spline->numberOfPoints - 1]));

    int whichInterval = binary_search(spline->numberOfPoints, spline->points, evaluationPoint);
    double pointsDifference = evaluationPoint - spline->points[whichInterval];
    double interpolantValue =
            spline->firstCoefficient[whichInterval] + 2 * pointsDifference * spline->secondCoefficient[whichInterval] +
            3 * pointsDifference * pointsDifference * spline->thirdCoefficient[whichInterval];

    return interpolantValue;
}

double evaluate_cubic_spline_integral(cubicSpline *spline, double evaluationPoints)
{
    //Binary search to find subinterval for evaluation point
    int whichInterval = binary_search((spline->numberOfPoints), spline->points, evaluationPoints);

    double integral = 0;
    double pointsDifference;

    for (int i = 0; i <= whichInterval; i++)
    {
        if (i < whichInterval)
        {
            pointsDifference = spline->points[i + 1] - spline->points[i];
        }
        else
        {
            pointsDifference = evaluationPoints - spline->points[i];
        }
        integral += spline->functionValues[i] * pointsDifference +
                    spline->firstCoefficient[i] * pointsDifference * pointsDifference / 2 +
                    spline->secondCoefficient[i] * pointsDifference * pointsDifference * pointsDifference / 3 +
                    spline->thirdCoefficient[i] * pointsDifference * pointsDifference * pointsDifference *
                    pointsDifference / 4;
    }
    return integral;
}

void free_cubic_spline(cubicSpline *spline)
{
    free(spline->points);
    free(spline->functionValues);
    free(spline->firstCoefficient);
    free(spline->secondCoefficient);
    free(spline);
}