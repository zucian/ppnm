#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include "utilities.h"
#include "integrateFunction.h"
#include "linearSpline.h"
#include "quadraticSpline.h"
#include "cubicSpline.h"

double gsl_cos(double x, void *params)
{
    return cos(x);
}

void cosine_integral(int numberOfPoints, double *xData, double *yData, double lowerLimit, double upperLimit,
                     double absoluteError,
                     double relativeError, size_t iterationLimit)
{
    double integralValue = linear_spline_integrate(numberOfPoints, xData, yData, upperLimit);

    //Gsl integrate settings
    gsl_function gslFuncCos;
    gslFuncCos.function = &gsl_cos;
    gslFuncCos.params = NULL;
    double integralValueGSL = integrate_function(lowerLimit, upperLimit, &gslFuncCos, absoluteError, relativeError,
                                                 iterationLimit);

    printf("Integral of interpolant from %g to %g: %g \n", lowerLimit, upperLimit, integralValue);
    printf("Integral of cos(x) from %g to %g using GSL: %g \n", lowerLimit, upperLimit, integralValueGSL);
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {    // Check that we have passed any arguments
        fprintf(stderr, "Error, no arguments were passed.\n"); // Else print to stderr
        exit(-1);
    }

    int numberOfPoints = 20;
    int numberOfSamples = (int) 1e3;

    //Insert data from cosine into arrays
    double *xData = malloc(numberOfPoints * sizeof(double));
    double *yData = malloc(numberOfPoints * sizeof(double));
    char *inputFile = argv[1];
    inputToArray(xData, yData, inputFile);

    //Define resolution of splines
    double resolution = fabs(xData[numberOfPoints - 1] - xData[0]) / numberOfSamples;

    double xDataTemp[3] = {xData[0],xData[1],xData[2]};

    printf("%g",xDataTemp[0]);
    printf("%g",xDataTemp[1]);
    printf("%g",xDataTemp[2]);

    free(xDataTemp);

    double xDataTemp[3] = {xData[1],xData[2],xData[3]};

    printf("%g",xDataTemp[0]);
    printf("%g",xDataTemp[1]);
    printf("%g",xDataTemp[2]);

    free(xDataTemp);

    printf("%g",xDataTemp[0]);
    printf("%g",xDataTemp[1]);
    printf("%g",xDataTemp[2]);



    //Build quadratic spline through each inner point with its neighbouring points
    /*
    for (int i = 2; i < xData[numberOfPoints - 1]; i += resolution)
    {
        for (int j = )
    }
    */

    return 0;
}