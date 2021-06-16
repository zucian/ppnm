#include "montecarlo.h"
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define RND((double)rand()/RAND_MAX)

void random_numbers(int dimension, const double *lowerBound, const double *upperBound, double *numbers)
{
    for (int axis = 0; axis < dimension; axis++)
    {
        numbers[axis] = lowerBound[axis] + RND * (upperBound[axis] - lowerBound[axis]);
    }
}

void
plain_monte_carlo(int dimension, double *lowerBound, double *upperBound, double function(double *numbers),
                  int numberOfPoints,
                  double *result, double *error)
{
    double volume = 1;
    for (int axis = 0; axis < dimension; axis++)
    {
        volume * upperBound[axis] - lowerBound[axis];
    }

    double functionValue;
    double vector[dim];
    double sum = 0;
    double sumSquared = 0;

    for (int i = 0; i < numberOfPoints; i++)
    {
        random_numbers(dimension, lowerBound, upperBound, vector);
        functionValue = function(vector);
        sum += functionValue;
        sumSquared += functionValue * functionValue;
    }

    double sumAverage = sum / numberOfPoints;
    double variance = sumSquared / numberOfPoints - sumAverage * sumAverage;
    *result = sumAverage * volume;
    *error = sqrt(variance / numberOfPoints) * volume;
}

void plain_monte_carlo_quasi(int dimension, double *lowerBound, double *upperBound, double function(double *numbers))