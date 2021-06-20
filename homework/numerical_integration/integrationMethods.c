#include<math.h>
#include<assert.h>
#include<stdio.h>
#include "integrationMethods.h"

double adapt_24(double function(double), double leftEndpoint, double rightEndpoint, double absoluteAccuracy,
                double relativeAccuracy, double secondFunctionValue, double thirdFunctionValue, int numberOfRecursions,
                double *integrationError)
{
    assert(numberOfRecursions < 1e6);

    double firstFunctionValue = function(leftEndpoint + (rightEndpoint - leftEndpoint) / 6);
    double fourthFunctionValue = function(leftEndpoint + 5 * (rightEndpoint - leftEndpoint) / 6);
    double higherOrderEvaluation =
            (2 * firstFunctionValue + secondFunctionValue + thirdFunctionValue + 2 * fourthFunctionValue) *
            (rightEndpoint - leftEndpoint) / 6;
    double lowerOrderEvaluation =
            (firstFunctionValue + secondFunctionValue + thirdFunctionValue + fourthFunctionValue) *
            (rightEndpoint - leftEndpoint) / 4;
    double tolerance = absoluteAccuracy + relativeAccuracy * fabs(higherOrderEvaluation);
    double error = fabs(higherOrderEvaluation - lowerOrderEvaluation);

    if (error < tolerance)
    {
        *integrationError += error;
        return higherOrderEvaluation;
    }
    else
    {
        double newHigherOrderEvaluationLeft = adapt_24(function, leftEndpoint, (leftEndpoint + rightEndpoint) / 2,
                                                       absoluteAccuracy / sqrt(2.), relativeAccuracy,
                                                       firstFunctionValue, secondFunctionValue, numberOfRecursions + 1,
                                                       integrationError);
        double newHigherOrderEvaluationRight = adapt_24(function, (leftEndpoint + rightEndpoint) / 2, rightEndpoint,
                                                        absoluteAccuracy / sqrt(2.), relativeAccuracy,
                                                        thirdFunctionValue, fourthFunctionValue, numberOfRecursions + 1,
                                                        integrationError);
        return newHigherOrderEvaluationLeft + newHigherOrderEvaluationRight;
    }

}

double adapt(double function(double), double leftEndpoint, double rightEndpoint, double absoluteAccuracy,
             double relativeAccuracy, double *integrationError)
{
    double secondFunctionValue = function(leftEndpoint + 2 * (rightEndpoint - leftEndpoint) / 6);
    double thirdFunctionValue = function(leftEndpoint + 4 * (rightEndpoint - leftEndpoint) / 6);
    int numberOfRecursion = 0;
    return adapt_24(function, leftEndpoint, rightEndpoint, absoluteAccuracy, relativeAccuracy, secondFunctionValue,
                    thirdFunctionValue, numberOfRecursion, integrationError);
}

double open_quad(double function(double), double leftEndpoint, double rightEndpoint, double absoluteAccuracy,
                 double relativeAccuracy, double *integrationError)
{
    double transformInterval(double x)
    {
        return 1.0 / 2 * (rightEndpoint - leftEndpoint) *
               function(1.0 / 2.0 * (leftEndpoint + rightEndpoint) + 1.0 / 2.0 * (rightEndpoint - leftEndpoint) * x);
    }
    double clenshawCurtisQuad(double angle)
    {
        return transformInterval(cos(angle)) * sin(angle);
    }
    return adapt(clenshawCurtisQuad, 0, M_PI, absoluteAccuracy, relativeAccuracy, integrationError);
}

double integrate(double function(double), double leftEndpoint, double rightEndpoint, double absoluteAccuracy,
                 double relativeAccuracy, double *integrationError)
{
    if (isinf(-leftEndpoint))
    {
        if (isinf(rightEndpoint))
        {
            double transformInterval(double t)
            {
                return function(t / (1 - t * t)) * (1 + t * t) / pow(1 - t * t, 2);
            }
            return adapt(transformInterval, -1, 1, absoluteAccuracy, relativeAccuracy, integrationError);
        }
        else
        {
            double transformInterval(double t)
            {
                return function(rightEndpoint + t / (1 + t)) / pow(1 + t, 2);
            }
            return adapt(transformInterval, -1, 0, absoluteAccuracy, relativeAccuracy, integrationError);
        }
    }
    else if (isinf(rightEndpoint))
    {
        double transformInterval(double t)
        {
            return function(leftEndpoint + t / (1 - t)) / pow(1 - t, 2);
        }
        return adapt(transformInterval, 0, 1, absoluteAccuracy, relativeAccuracy, integrationError);
    }
    else
    {
        return adapt(function, leftEndpoint, rightEndpoint, absoluteAccuracy, relativeAccuracy, integrationError);
    }
}