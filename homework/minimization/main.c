#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "minimization.h"
#include "inputArray.h"

double energy;
double boundEnergy;
double maxPoint;
double boundMaxPoint;
double *error;
int numberOfDataPoints;
double *energyBW;
double *crossSection;

double simplex_test_function(double *values)
{
    double x = values[0];
    double y = values[0];

}

double rosenbrock_valley(gsl_vector *values)//Minimum at (a,a^2) = (1,1)
{
    double x = gsl_vector_get(values, 0);
    double y = gsl_vector_get(values, 1);
    double functionValue = (1 - x) * (1 - x) + 100 * (y - x * x) * (y - x * x);
    return functionValue;
}

double himmelblau(gsl_vector *values)
{
    double x = gsl_vector_get(values, 0);
    double y = gsl_vector_get(values, 1);
    double functionValue = (x * x + y - 11) * (x * x + y - 11) + (x + y * y - 7) * (x + y * y - 7);
    return functionValue;
}

double breit_Wigner(double mass, double width, double scaleFactor, double energy)
{
    return scaleFactor / ((energy - mass) * (energy - mass) + (width * width) / 4);
}

double deviation_function(gsl_vector *values)
{
    double mass = gsl_vector_get(values, 0);
    double width = gsl_vector_get(values, 1);
    double scaleFactor = gsl_vector_get(values, 2);
    double sum = 0;

    for (int i = 0; i < numberOfDataPoints; i++)
    {
        sum += pow(breit_Wigner(mass, width, scaleFactor, energyBW[i]) - crossSection[i], 2) / pow(error[i], 2);
    }
}

int main(int argc, char *argv[])
{
    //PART A
    printf("A: quasi-Newton minimization with numerical gradient, back-tracking linesearch, rank-1 update\n\n");
    int dimension = 2;
    double tolerance = 1e-5;
    gsl_vector *minimum = gsl_vector_alloc(dimension);

    //ROSENBROCK
    printf("Testing minimization routine on Rosenbrock's valley function\n");

    double minimumX = 1; //x-value of minimum for Rosenbrock
    double minimumY = 1; //y-value of minimum for Rosenbrock
    double initialXValue = minimumX - 1;
    double initialYValue = minimumY - 1;

    gsl_vector_set(minimum, 0, initialXValue);
    gsl_vector_set(minimum, 0, initialYValue);

    printf("Initial value (x,y): (%g,%g)\n", initialXValue, initialYValue);
    quasi_newton_method(rosenbrock_valley, minimum, tolerance);
    printf("Found minimum (x,y): (%g,%g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    printf("Actual global minimum is (1,1)\n\n");

    //Himmelblaus
    printf("Testing minimization routine on Himmelblau's function\n");

    minimumX = 3;
    minimumY = 2;
    initialXValue = minimumX - 0.2;
    initialYValue = minimumY - 0.2;

    gsl_vector_set(minimum, 0, initialXValue);
    gsl_vector_set(minimum, 0, initialYValue);

    printf("Initial value (x,y): (%g,%g)\n", initialXValue, initialYValue);
    quasi_newton_method(himmelblau, minimum, tolerance);
    printf("Found minimum (x,y): (%g,%g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    printf("Actual minima are four points (x,y): (%g, %g), (-3.78, -3.28), (-2.80, 3.13), (3.58, -1.84)\n\n", minimumX,
           minimumY);

    //PART B

    printf("B: Higgs discovery\n\n");

    double measuredMass = 125.3;
    numberOfDataPoints = 30;
    energyBW = (double *) malloc(numberOfDataPoints * sizeof(double));
    crossSection = (double *) malloc(numberOfDataPoints * sizeof(double));
    error = (double *) malloc(numberOfDataPoints * sizeof(double));

    inputArray(energyBW, crossSection, error, "higgsData.txt");
    printf("energy E[GeV], cross section σ(E), error δσ \n");

    for (int i = 0; i < numberOfDataPoints; i++)
    {
        printf("%g \t %g \t %g\n\n", energyBW[i], crossSection[i], error[i]);
    }

    dimension = 3;
    gsl_vector *minimumBW = gsl_vector_alloc(dimension);
    gsl_vector_set(minimumBW, 0, measuredMass + 1);
    gsl_vector_set(minimumBW, 1, 2.7);
    gsl_vector_set(minimumBW, 2, 8);

    quasi_newton_method(deviation_function, minimumBW, tolerance);
    printf("Initial value (m, Γ, A): (%g, %g, %g) \n", measuredMass + 1, 2.7, 8.0);
    printf("Found Minima (m, Γ, A): (%g, %g, %g) \n\n", gsl_vector_get(minimumBW, 0), gsl_vector_get(minimumBW, 1),
           gsl_vector_get(minimumBW, 2));

    FILE *outputFilestream = fopen("higgsFit.txt", "w");
    for (int i = 0; i < numberOfDataPoints; i++)
    {
        double functionValue = breit_Wigner(gsl_vector_get(minimumBW, 0), gsl_vector_get(minimumBW, 1),
                                            gsl_vector_get(minimumBW, 2), energyBW[i]);
        fprintf(outputFilestream, "%g \t %g \t %g \t %g \n", energyBW[i], crossSection[i], error[i], functionValue);
    }
    fclose(outputFilestream);

    printf("Fit can be seen in file higgsFit.png\n");


    //PART C
    printf("C: Implement downhill simplex method");



    return 0;
}