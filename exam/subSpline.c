#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "utilities.h"
#include "subSpline.h"
#include "quadraticSpline.h"

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

    double deltaX;
    double deltaY;
    double deltaP;
    //Compute c_i and d_i
    for (int i = 0; i < numberOfPoints - 1; i++)
    {
        //Defining Δx, Δy, Δp as difference between neighbour points
        deltaX = spline->points[i + 1] - spline->points[i];
        deltaY = spline->functionValues[i + 1] - spline->functionValues[i];
        deltaP = spline->derivativeValues[i + 1] - spline->derivativeValues[i];

        //c_i (second coefficient) and d_i (third coefficient)
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

    double interpolantValue = (spline->functionValues[whichInterval]) + (spline->firstCoefficient[whichInterval]) *
                                                                        (evaluationPoint -
                                                                         spline->points[whichInterval]) +
                              (spline->secondCoefficient[whichInterval]) *
                              (evaluationPoint - spline->points[whichInterval]) *
                              (evaluationPoint - spline->points[whichInterval]) +
                              (spline->thirdCoefficient[whichInterval]) *
                              (evaluationPoint - spline->points[whichInterval]) *
                              (evaluationPoint - spline->points[whichInterval]) *
                              (evaluationPoint - spline->points[whichInterval]);

    return interpolantValue;
}

void estimate_derivative(int numberOfPoints, double *xData, double *yData, double *pData)
{
    //The function will fill pData with estimates of first derivative using quadratic splines

    //Amount of points for quadratic spline to estimate derivative
    int temporaryNumberOfPoints = 3;

    double xDataTemporary[3]; //Array that holds points x_{i-1}, x_i, x_{i+1}
    double yDataTemporary[3]; //Array that holds points y_{i-1}, y_i, y_{i+1}

    int counter = 1; //Counter that checks whether we are at first or last point

    for (int i = 1; i < numberOfPoints - 1; i++)
    {
        //Fill subarrays with data
        xDataTemporary[0] = xData[i - 1];
        xDataTemporary[1] = xData[i];
        xDataTemporary[2] = xData[i + 1];
        yDataTemporary[0] = yData[i - 1];
        yDataTemporary[1] = yData[i];
        yDataTemporary[2] = yData[i + 1];

        //Create quadratic spline using temporary data
        quadSpline *temporaryQuadraticSpline = initialize_quadratic_spline(temporaryNumberOfPoints, xDataTemporary,
                                                                           yDataTemporary);

        //If this is first point, use polynomial to estimate p_1
        if (counter == 1)
        {
            double temporaryInterpolantDerivativeStart = evaluate_quadratic_spline_derivative(temporaryQuadraticSpline,
                                                                                              xDataTemporary[0]);
            pData[i - 1] = temporaryInterpolantDerivativeStart;
        }

        //Evaluate p_i
        double temporaryInterpolantDerivative = evaluate_quadratic_spline_derivative(temporaryQuadraticSpline,
                                                                                     xDataTemporary[1]);
        pData[i] = temporaryInterpolantDerivative;

        counter++;

        //If this is last point, use polynomial to estimate p_n
        if (counter == numberOfPoints - 1)
        {
            double temporaryInterpolantDerivativeEnd = evaluate_quadratic_spline_derivative(temporaryQuadraticSpline,
                                                                                            xDataTemporary[2]);
            pData[i + 1] = temporaryInterpolantDerivativeEnd;
        }

    }
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
