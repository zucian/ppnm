#include <assert.h>
#include "linearSpline.h"
#include "utilities.h"

double linear_spline(int numberOfPoints, double *points, double *functionValues, double evaluationPoints)
{
    //Ensure proper points
    assert((numberOfPoints) > 1 && (evaluationPoints >= points[0]) && (evaluationPoints <= points[numberOfPoints]));

    //Binary search to find subinterval for evaluation point
    int whichInterval = binary_search(numberOfPoints, points, evaluationPoints);
    double functionValueDifference = (functionValues[whichInterval + 1] - functionValues[whichInterval]);
    double pointsDifference = (points[whichInterval + 1] - points[whichInterval]);
    double slope = functionValueDifference / pointsDifference;

    assert(pointsDifference > 0);
    double interpolationValue = functionValues[whichInterval] + slope * (evaluationPoints - points[whichInterval]);
    return interpolationValue;
}

double linear_spline_integrate(int numberOfPoints, double *points, double *functionValues, double evaluationPoint)
{
    //Binary search to find subinterval for evaluation point
    int whichInterval = binary_search(numberOfPoints, points, evaluationPoint);

    double integral = 0;
    double functionValueDifference;
    double pointsDifference;
    double slope;

    for (int i = 0; i <= whichInterval; i++)
    {
        functionValueDifference = functionValues[i + 1] - functionValues[i];
        pointsDifference = points[i + 1] - points[i];
        slope = functionValueDifference / pointsDifference;

        if (i >= whichInterval)
        {
            pointsDifference = evaluationPoint - points[i];
        }
        //analytical solution of integral
        integral += functionValues[i] * pointsDifference + slope * pointsDifference * pointsDifference / 2;
    }
    return integral;
}
