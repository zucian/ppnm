#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "utilities.h"
#include "subSpline.h"

subSpline *initialize_sub_spline(int numberOfPoints, double *xData, double *yData, double *pData)
{
    //Allocate memory for subspline
    subSpline *spline = (subSpline *) malloc(sizeof(subSpline));
    spline->points = (double *) malloc(numberOfPoints * sizeof(double));
    spline->functionValues = (double *) malloc(numberOfPoints * sizeof(double));
    spline->derivativeValues = (double *) malloc(numberOfPoints * sizeof(double));
    spline->firstCoefficient = (double *) malloc(numberOfPoints * sizeof(double));
    spline->secondCoefficient = (double *) malloc(numberOfPoints * sizeof(double));
    spline->thirdCoefficient = (double *) malloc(numberOfPoints * sizeof(double));
    spline->numberOfPoints = numberOfPoints;

    //Fill spline values from input data
    for (int i = 0; i < numberOfPoints; i++)
    {
        spline->points[i] = xData[i];
        spline->functionValues[i] = yData[i];
        spline->derivativeValues[i] = pData[i];
    }

    //Set b_i (first coefficient) equal to p_i
    for (int i = 0; i < numberOfPoints; i++)
    {
        spline->firstCoefficient[i] = spline->derivativeValues[i];
    }

    //Compute c_i and d_i
    for (int i = 0; i < numberOfPoints - 1; i++)
    {
        //Defining Δx, Δy, Δp as difference between neighbour points
        double deltaX = spline->points[i + 1] - spline->points[i];
        double deltaY = spline->functionValues[i + 1] - spline->functionValues[i];
        double deltaP = spline->derivativeValues[i + 1] - spline->derivativeValues[i];

        //c_i (second coefficient) and d_i (thrid coefficient)
        spline->secondCoefficient[i] =
                3 * deltaY / (deltaX * deltaX) - 3 * spline->derivativeValues[i] / deltaX - deltaP / deltaX;
        spline->thirdCoefficient[i] = deltaP / (deltaX * deltaX) - 2 * deltaY / (deltaX * deltaX * deltaX) +
                                      2 * spline->derivativeValues[i] / (deltaX * deltaX);
    }

    return spline;
}

double evaluate_sub_spline(subSpline *spline, double evaluationPoint)
{
    //Make sure we evaluate spline inside interval for data points
    assert((evaluationPoint >= (spline->points[0])) &&
           (evaluationPoint <= (spline->points[spline->numberOfPoints - 1])));

    //Binary search to find subinterval for evaluation point
    int whichInterval = binarySearch(spline->numberOfPoints, spline->points, evaluationPoint);

    double pointsDifference = evaluationPoint - (spline->points[whichInterval]);
    double thirdDifference = pointsDifference * (spline->thirdCoefficient[whichInterval]);
    double secondDifference = pointsDifference * (spline->secondCoefficient[whichInterval]);
    double firstDifference = pointsDifference * (spline->thirdCoefficient[whichInterval]);

    double interpolantValue = (spline->functionValues[whichInterval]) + firstDifference;

    return interpolantValue;
}

void free_sub_spline(subSpline *spline)
{
    free(spline->points);
    free(spline->functionValues);
    free(spline->derivativeValues);
    free(spline->firstCoefficient);
    free(spline->secondCoefficient);
    free(spline->thirdCoefficient);
}
