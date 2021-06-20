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
    double secondDerivative = (-2) * (1.0 / variable + energy) * currentFunctionValue;
    gsl_vector_set(functionDerivative, 0, firstDerivative);
    gsl_vector_set(functionDerivative, 1, secondDerivative);
}

void schrodinger_bound(double variable, gsl_vector *functionValues, gsl_vector *functionDerivative)
{
    double currentFunctionValue = gsl_vector_get(functionValues, 0);
    double firstDerivative = gsl_vector_get(functionValues, 1);
    double secondDerivative = (-2) * (1.0 / variable + boundEnergy) * currentFunctionValue;
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
    //PART A
    printf("A: Newton's method with numerical Jacobian and back-tracking linesearch \n\n");

    printf("Testing routine on function (x-4)^2+(y-9)^2+1\n\n");

    int numberOfDimensions = 2;
    double tolerance = 1e-5;

    //Set initial values for guess
    gsl_vector *minimum = gsl_vector_alloc(numberOfDimensions);
    gsl_vector_set(minimum, 0, 3);
    gsl_vector_set(minimum, 1, 7);

    printf("Initial guess (x,y): (%g,%g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    newton_method(test_function, minimum, tolerance);
    printf("Found minimum (x,y): (%g,%g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    printf("Actual minimum is (4,9) \n");
    printf("Tolerance used was %g\n\n", tolerance);

    printf("Testing routine on Rosenbrock valley function \n\n");

    //Set initial values for guess
    gsl_vector_set(minimum, 0, 0.5);
    gsl_vector_set(minimum, 1, 0.5);

    printf("Initial guess (x,y): (%g,%g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    newton_method(rosenbrock_gradient, minimum, tolerance);
    printf("Found minimum (x,y): (%g,%g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    printf("Actual minimum is (1,1) \n");
    printf("Tolerance used was %g\n\n", tolerance);

    //PART B
    printf("B: Bound state of hydrogen atom \n\n");

    gsl_vector *minimumHydrogen = gsl_vector_alloc(1);
    gsl_vector *minimumHydrogenBound = gsl_vector_alloc(1);
    gsl_vector_set(minimumHydrogen, 0, -3);
    gsl_vector_set(minimumHydrogenBound, 0, -1);
    maxPoint = 8.0;
    boundMaxPoint = 0.5;

    //Finding energies
    newton_method(wavefunction, minimumHydrogen, tolerance);
    newton_method(wavefunction_bound, minimumHydrogenBound, tolerance);

    printf("Unbound energy found using root-finding: %g\n", energy);

    //Settings for ODE
    int dimension = 2;
    gsl_vector *functionValueLeft = gsl_vector_calloc(dimension);
    gsl_vector *functionValueRight = gsl_vector_alloc(dimension);
    double leftEndpoint = 1e-3;
    double rightEndpoint = maxPoint;
    double absoluteAccuracy = 1e-3;
    double relativeAccuracy = 1e-3;
    double step = (rightEndpoint - leftEndpoint) / 10;

    gsl_vector_set(functionValueLeft, 0, (leftEndpoint - leftEndpoint * leftEndpoint));
    gsl_vector_set(functionValueLeft, 1, (1 - 2. * leftEndpoint));

    printf("Solving ODE with found energy: \n");
    FILE *ODEOutput = fopen(argv[1], "w");
    rk_driver(schrodinger, leftEndpoint, functionValueLeft, rightEndpoint, functionValueRight, step, absoluteAccuracy,
              relativeAccuracy, ODEOutput);
    fclose(ODEOutput);
    printf("Result plotted in hydrogenPlot.png\n\n");

    //PART C
    printf("C: Better boundary condition for hydrogen atom problem\n\n");

    printf("Bound energy found using root-finding: %g \n\n", boundEnergy); //Done in part B

    printf("Investigating convergence of energy minimum as function of r_max \n\n");

    double exactEnergy = -1.0 / 2.0;
    FILE *convergenceData = fopen("convergenceData.txt", "w");
    for (double i = 0.1; i <= 8.0; i += 0.1)
    {
        maxPoint = i;
        boundMaxPoint = i;

        //Reset minimumHydrogen
        gsl_vector_set(minimumHydrogen, 0, -3);
        gsl_vector_set(minimumHydrogenBound, 0, -1);

        newton_method(wavefunction, minimumHydrogen, tolerance);
        newton_method(wavefunction_bound, minimumHydrogenBound, tolerance);

        fprintf(convergenceData, "%g \t %g \t %g \n", maxPoint, fabs(energy-exactEnergy), fabs(boundEnergy-exactEnergy));
    }
    fclose(convergenceData);
    printf("Convergence results can be seen in convergencePlot.png \n");

    //Free memory
    gsl_vector_free(minimum);
    gsl_vector_free(minimumHydrogen);
    gsl_vector_free(minimumHydrogenBound);

    return 0;
}