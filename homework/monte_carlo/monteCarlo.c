#include "monteCarlo.h"
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define fracl(x) ((x) -  floorl(x))
#define PI 3.1415926535897932384626433832795028841971693993751L
#define RND ((double)rand()/RAND_MAX)

void random_numbers(int dimension, const double *lowerBound, const double *upperBound, double *numbers)
{
    for (int axis = 0; axis < dimension; axis++)
    {
        numbers[axis] = lowerBound[axis] + RND * (upperBound[axis] - lowerBound[axis]);
    }
}

double van_der_corput_sequence(int i, int base)
{
    double corputNumber = 0;
    double coprimeBase = (double) 1 / base;
    while (i > 0)
    {
        corputNumber += (i % base) * coprimeBase;
        i /= base;
        coprimeBase /= base;
    }
    return corputNumber;
}

double halton_sequence_first(int i, int dimension, double *vector)
{
    int base[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43};
    int maximumDimension = sizeof(base) / sizeof(int);

    assert(dimension <= maximumDimension);

    for (int axis = 0; axis < dimension; axis++)
    {
        vector[axis] = van_der_corput_sequence(i + 1, base[axis]);
    }
}

double halton_sequence_second(int i, int dimension, double *vector)
{
    int base[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
    int maximumDimension = sizeof(base) / sizeof(int);

    assert(dimension <= maximumDimension);

    for (int axis = 0; axis < dimension; axis++)
    {
        vector[axis] = van_der_corput_sequence(i + 1, base[axis]);
    }
}

void random_number_halton_corput_first(int i, int dimension, const double *lowerBound, const double *upperBound,
                                       double *vector)
{
    halton_sequence_first(i, dimension, vector);
    for (int axis = 0; axis < dimension; axis++)
    {
        vector[axis] = lowerBound[axis] + (upperBound[axis] - lowerBound[axis]) * vector[axis];
    }
}

void random_number_halton_corput_second(int i, int dimension, const double *lowerBound, const double *upperBound,
                                        double *vector)
{
    halton_sequence_second(i, dimension, vector);
    for (int axis = 0; axis < dimension; axis++)
    {
        vector[axis] = lowerBound[axis] + (upperBound[axis] - lowerBound[axis]) * vector[axis];
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
        volume *= upperBound[axis] - lowerBound[axis];
    }

    double functionValue;
    double vector[dimension];
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

void plain_monte_carlo_quasi(int dimension, double *lowerBound, double *upperBound, double function(double *numbers),
                             int numberOfPoints, double *result, double *error)
{
    double volume = 1;
    for (int axis = 0; axis < dimension; axis++)
    {
        volume *= upperBound[axis] - lowerBound[axis];
    }

    double functionValueFirst;
    double functionValueSecond;
    double vectorFirst[dimension];
    double vectorSecond[dimension];
    double sumFirst;
    double sumSecond;

    for (int i = 0; i < floor(numberOfPoints / 2); i++)
    {
        random_number_halton_corput_first(i, dimension, lowerBound, upperBound, vectorFirst);
        random_number_halton_corput_second(i, dimension, lowerBound, upperBound, vectorSecond);
        functionValueFirst = function(vectorFirst);
        functionValueSecond = function(vectorSecond);

        if (!isinf(functionValueFirst) && !isinf(functionValueSecond))
        {
            sumFirst += functionValueFirst;
            sumSecond += functionValueSecond;
        }
    }

    double average = (sumFirst + sumSecond) / numberOfPoints;

    *result = volume * average;
    *error = volume * fabs(sumFirst - sumSecond) / numberOfPoints;
}

double plain_monte_carlo_stratified_sampling(int dimension, double function(double *vector), double *lowerBound,
                                             double *upperBound, double absoluteAccuracy, double relativeAccuracy,
                                             int numberOfRecalls, double meanRecalls)
{
    int numberOfPoints = 16 * dimension;
    double volume = 1;

    for (int axis = 0; axis < dimension; axis++)
    {
        volume *= upperBound[axis] - lowerBound[axis];
    }

    double average = 0;
    double averageLeft[dimension];
    double averageRight[dimension];
    double vector[dimension];
    int numberLeft[dimension];
    int numberRight[dimension];

    for (int axis = 0; axis < dimension; axis++)
    {
        averageLeft[axis] = 0;
        averageRight[axis] = 0;
        numberLeft[axis] = 0;
        numberRight[axis] = 0;
    }

    for (int i = 0; i < numberOfPoints; i++)
    {
        random_numbers(dimension, lowerBound, upperBound, vector);
        double functionValue = function(vector);
        average += functionValue;

        for (int axis = 0; axis < dimension; axis++)
        {
            double middlePoint = (lowerBound[axis] + upperBound[axis]) / 2;

            if (vector[axis] > middlePoint)
            {
                averageRight[axis] += functionValue;
                numberRight[axis]++;
            }
            else
            {
                averageLeft[axis] += functionValue;
                numberLeft[axis]++;
            }
        }
    }

    average /= numberOfPoints;

    for (int axis = 0; axis < dimension; axis++)
    {
        averageLeft[axis] /= numberLeft[axis];
        averageRight[axis] /= numberRight[axis];
    }

    int div = 0;
    double maxVariance = 0;

    for (int axis = 0; axis < dimension; axis++)
    {
        double variance = fabs(averageRight[axis] - averageLeft[axis]);
        if (variance > maxVariance)
        {
            maxVariance = variance;
            div = axis;
        }
    }

    //Final results
    double result =
            (average * numberOfPoints + meanRecalls * numberOfRecalls) / (numberOfPoints + numberOfRecalls) * volume;
    double error = volume * fabs(meanRecalls - average);
    double tolerance = absoluteAccuracy + relativeAccuracy * fabs(result);


    if (error < tolerance)
    {
        return result;
    }

    double lowerBoundSecond[dimension];
    double upperBoundSecond[dimension];

    for (int axis = 0; axis < dimension; axis++)
    {
        upperBoundSecond[axis] = upperBound[axis];
    }

    lowerBoundSecond[div] = (lowerBound[div] + upperBound[div]) / 2;
    upperBoundSecond[div] = (lowerBound[div] + upperBound[div]) / 2;

    //Call recursively
    double resultLeft = plain_monte_carlo_stratified_sampling(dimension, function, lowerBound, upperBoundSecond,
                                                              absoluteAccuracy / sqrt(2), relativeAccuracy,
                                                              numberLeft[div], averageLeft[div]);
    double resultRight = plain_monte_carlo_stratified_sampling(dimension, function, lowerBoundSecond, upperBound,
                                                               absoluteAccuracy / sqrt(2), relativeAccuracy,
                                                               numberRight[div], averageRight[div]);

    return resultLeft + resultRight;
}

void lattice(int d, double *x)
{
    static int dimension = 0;
    static int n = -1;
    static long double *alpha;

    if (d < 0)
    {
        dimension = -d;
        n = 0;
        alpha = (long double *) realloc(alpha, dimension * sizeof(long double));

        for (int i = 0; i < dimension; i++)
        {
            alpha[i] = fracl(sqrtl(PI + 1));
        }
    }
    else if (d > 0)
    {
        n++;
        assert(d == dimension && n > 0);
        for (int i = 0; i < dimension; i++)
        {
            x[i] = fracl(n * alpha[i]);
        }
    }
    else if (alpha != NULL)
    {
        free(alpha);
    }
}