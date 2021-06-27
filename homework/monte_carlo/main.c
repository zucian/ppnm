#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "monteCarlo.h"

double debug_function(double *x) //sqrt(x) function for testing monte carlo integration
{
    return sqrt(x[0]);
}

double dmitri_function(double *x) //Function supplied in homework description
{
    return 1 / (M_PI * M_PI * M_PI * (1 - cos(x[0]) * cos(x[1]) * cos(x[2])));
}

int main(int argc, char *argv[])
{
    int numberOfPoints = (int) 1e6;

    // PART A
    printf("\nA: Plain Monte Carlo Integration \n\n");

    //Square root test function
    printf("Testing the plain Monte Carlo routine on √(x) from 0 to 2 \n");
    double exactValueTestFunction = 1.8856180831641;
    printf("Exact value of integral is %g\n", exactValueTestFunction);

    const int dimensionDebug = 1;
    double *lowerBoundDebug = calloc(dimensionDebug, sizeof(double));
    double *upperBoundDebug = malloc(sizeof(double) * dimensionDebug);
    lowerBoundDebug[0] = 0;
    upperBoundDebug[0] = 2;
    double resultDebug = 0;
    double errorDebug = 0;

    plain_monte_carlo(dimensionDebug, lowerBoundDebug, upperBoundDebug, debug_function, numberOfPoints, &resultDebug,
                      &errorDebug);
    printf("Numerical estimate using Monte Carlo: \n");
    printf("Estimate: %g\n", resultDebug);
    printf("Error: %g\n", fabs(resultDebug - exactValueTestFunction));
    printf("Error estimate from monte carlo: %g\n\n", errorDebug);

    //The required function
    printf("Testing the plain Monte Carlo routine on ∫_0^π dx/π ∫_0^π dy/π ∫_0^π  dz/π [1-cos(x)cos(y)cos(z)]^{-1} = Γ(1/4)4/(4π3)\n");
    double exactValueDmitriFunction = 1.3932039296856768591842462603255;
    printf("Exact value of integral is %g\n", exactValueDmitriFunction);

    const int dimensionDmitri = 3;
    double *lowerBoundDmitri = calloc(dimensionDmitri, sizeof(double));
    double *upperBoundDmitri = malloc(sizeof(double) * dimensionDmitri);
    for (int axis = 0; axis < dimensionDmitri; axis++)
    {
        lowerBoundDmitri[axis] = 0;
        upperBoundDmitri[axis] = M_PI;
    }

    double resultDmitri = 0;
    double errorDmitri = 0;

    plain_monte_carlo(dimensionDmitri, lowerBoundDmitri, upperBoundDmitri, dmitri_function, numberOfPoints,
                      &resultDmitri,
                      &errorDmitri);

    printf("Numerical estimate using Monte Carlo: \n");
    printf("Estimate: %g\n", resultDmitri);
    printf("Error: %g\n", fabs(resultDmitri - exactValueDmitriFunction));
    printf("Error estimate from monte carlo: %g\n\n", errorDmitri);

    //PART B
    printf("B: Quasi Monte Carlo Integration\n\n");

    double resultDmitriQuasi = 0;
    double errorDmitriQuasi = 0;

    plain_monte_carlo_quasi(dimensionDmitri, lowerBoundDmitri, upperBoundDmitri, dmitri_function, numberOfPoints,
                            &resultDmitriQuasi,
                            &errorDmitriQuasi);
    printf("Testing the quasi Monte Carlo routine on ∫_0^π dx/π ∫_0^π dy/π ∫_0^π  dz/π [1-cos(x)cos(y)cos(z)]^{-1} = Γ(1/4)4/(4π3)\n");
    printf("Numerical estimate using quasi Monte Carlo: \n");
    printf("Estimate: %g\n", resultDmitriQuasi);
    printf("Error: %g\n", fabs(resultDmitriQuasi - exactValueDmitriFunction));
    printf("Error estimate from monte carlo quasi: %g\n\n", errorDmitriQuasi);

    printf("Test error scaling of quasi- vs pseudo-random \n");

    char filename[] = "errorScaling.txt";
    printf("Result in errorScaling.png \n\n");
    int numberOfReps = (int) 500;
    int stepSize = (int) 100;

    FILE *outputFileStream = fopen(filename, "w");
    for (int rep = 1; rep <= numberOfReps; rep++)
    {
        int repNumberOfPoints = rep * stepSize;
        plain_monte_carlo(dimensionDmitri, lowerBoundDmitri, upperBoundDmitri, dmitri_function, repNumberOfPoints,
                          &resultDmitri,
                          &errorDmitri);
        plain_monte_carlo_quasi(dimensionDmitri, lowerBoundDmitri, upperBoundDmitri, dmitri_function, repNumberOfPoints,
                                &resultDmitriQuasi,
                                &errorDmitriQuasi);

        fprintf(outputFileStream, "%i \t %g \t %g \n", repNumberOfPoints, errorDmitri, errorDmitriQuasi);
    }
    fclose(outputFileStream);
    /*
    //PART C
    printf("C: Recursive stratified sampling\n\n");

    double absoluteAccuracy = 1e-3;
    double relativeAccuracy = 1e-3;
    double errorDmitriStratified = 0;

    int numberOfRecalls = 0;
    double meanRecalls = 0;
    double resultDmitriStratified = plain_monte_carlo_stratified_sampling(dimensionDmitri, dmitri_function,
                                                                          lowerBoundDmitri, upperBoundDmitri,
                                                                          absoluteAccuracy, relativeAccuracy,
                                                                          numberOfRecalls, meanRecalls);

    printf("Numerical estimate using Monte Carlo Stratified: \n");
    printf("Estimate: %g\n", resultDmitriStratified);
    printf("Error: %g\n", fabs(resultDmitriStratified - exactValueDmitriFunction));
    */
    return 0;
}
