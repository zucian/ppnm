#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "GS_utilities.h"
#include "rootFindingMethods.h"
#include "rungeKutta.h"

double energy;
double boundEnergy;
double maxPoint;
double boundMaxPoint;

void rosenbrock_gradient(gsl_vector *values, gsl_vector *functionValues)
{
    double x = gsl_vector_get(values, 0);
    double y = gsl_vector_get(values, 1);
    gsl_vector_set(functionValues, 0, 2 * (-1) * (1 - x) + (-2 * x) * 2 * 100 * (y - x * x));
    gsl_vector_set(functionValues, 1, 2 * 100 * (y - x * x));
}

void test_function(gsl_vector *values, gsl_vector *functionValues)
{
    double x = gsl_vector_get(values, 0);
    double y = gsl_vector_get(values, 1);
    double a = 4;
    double b = 9;
    gsl_vector_set(functionValues, 0, 2 * x * (x - a));
    gsl_vector_set(functionValues, 1, 2 * y * (y - b));
}

void schrodinger(double variable, gsl_vector *functionValues, gsl_vector *functionDerivative)
{
    double currentFunctionValue = gsl_vector_get(functionValues, 0);
    double firstDerivative = gsl_vector_get(functionValues, 1);
    double secondDerivative = (-2) * (1.0 / var + energy) * currentFunctionValue;
    gsl_vector_set(functionDerivative, 0, firstDerivative);
    gsl_vector_set(functionDerivative, 1, secondDerivative);
}

void schrodinger_bound(double variable, gsl_vector *functionValues, gsl_vector *functionDerivative)
{
    double currentFunctionValue = gsl_vector_get(functionValues, 0);
    double firstDerivative = gsl_vector_get(functionValues, 1);
    double secondDerivative = (-2) * (1.0 / var + boundEnergy) * currentFunctionValue;
    gsl_vector_set(functionDerivative, 0, firstDerivative);
    gsl_vector_set(functionDerivative, 1, secondDerivative);
}

void wavefunction(gsl_vector *values, gsl_vector *functionValues)
{
    int dimension = 2;
    gsl_vector *functionValueLeft = gsl_vector_alloc(dimension);
    gsl_vector *functionValueRight = gsl_vector_alloc(dimension);

    double leftEndpoint = 1e-3; //Cant choose zero because of divergence at zero, use small value instead
    double rightEndpoint = maxPoint;
    double absoluteAccuracy = 1e-3;
    double relativeAccuracy = 1e-3;
    double step = (rightEndpoint - leftEndpoint) / 10;
    energy = gsl_vector_get(values, 0);

    gsl_vector_set(functionValueLeft, 0, (leftEndpoint - leftEndpoint * leftEndpoint));
    gsl_vector_set(functionValueLeft, 1, (1 - 2 * leftEndpoint));

    rk_driver(schrodinger, leftEndpoint, functionValueLeft, rightEndpoint, functionValueRight, step, absoluteAccuracy,
              relativeAccuracy, NULL);
    gsl_vector_set(functionValues, 0, gsl_vector_get(functionValueRight, 0));
}

void wavefunction_bound(gsl_vector *values, gsl_vector *functionValues)
{
    int dimension = 2;
    gsl_vector *functionValueLeft = gsl_vector_alloc(dimension);
    gsl_vector *functionValueRight = gsl_vector_alloc(dimension);

    double leftEndpoint = 1e-5; //Cant choose zero because of divergence at zero, use small value instead
    double rightEndpoint = boundMaxPoint;
    double absoluteAccuracy = 1e-5;
    double relativeAccuracy = 1e-5;
    double step = (rightEndpoint - leftEndpoint) / 10;
    boundEnergy = gsl_vector_get(values, 0);

    gsl_vector_set(functionValueLeft, 0, (leftEndpoint - leftEndpoint * leftEndpoint));
    gsl_vector_set(functionValueLeft, 1, (1 - 2 * leftEndpoint));

    rk_driver(schrodinger_bound, leftEndpoint, functionValueLeft, rightEndpoint, functionValueRight, step,
              absoluteAccuracy,
              relativeAccuracy, NULL);
    gsl_vector_set(functionValues, 0, gsl_vector_get(functionValueRight, 0) -
                                      rightEndpoint * exp(-sqrt((-2) * boundEnergy) * rightEndpoint));
}

int main(int argc, char *argv[])
{
    return 0;
}