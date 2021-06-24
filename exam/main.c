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

    //Open files for data
    char *inputFile = argv[1];
    FILE *linearOutput = fopen(argv[2], "w");
    FILE *quadraticOutput = fopen(argv[3], "w");
    FILE *cubicOutput = fopen(argv[4], "w");

    double *xData = malloc(numberOfPoints * sizeof(double));
    double *yData = malloc(numberOfPoints * sizeof(double));

    inputToArray(xData, yData, inputFile);

    //Build quadratic spline through each inner point with its neighbouring points
    int temporaryNumberOfPoints = 3;

    double xDataTemporary[3]; //Array that holds points x_{i-1}, x_i, x_{i+1}
    double yDataTemporary[3]; //Array that holds points y_{i-1}, y_i, y_{i+1}

    xDataTemporary[0] = 0.0;
    xDataTemporary[1] = 0.0;
    xDataTemporary[2] = 0.0;
    yDataTemporary[0] = 0.0;
    yDataTemporary[1] = 0.0;
    yDataTemporary[2] = 0.0;

    quadSpline *temporaryQuadraticSpline = initialize_quadratic_spline(temporaryNumberOfPoints, xDataTemporary,
                                                                       yDataTemporary);

    for (int i = 1; i < numberOfPoints-1; i++)
    {
        //Fill subarrays with data
        xDataTemporary[0] = xData[i - 1];
        xDataTemporary[1] = xData[i];
        xDataTemporary[2] = xData[i + 1];
        yDataTemporary[0] = yData[i - 1];
        yDataTemporary[1] = yData[i];
        yDataTemporary[2] = yData[i + 1];

        //Create quadratic spline settings
        temporaryQuadraticSpline = initialize_quadratic_spline(temporaryNumberOfPoints, xDataTemporary, yDataTemporary);

        double temporaryInterpolantDerivative = evaluate_quadratic_spline_derivative(temporaryQuadraticSpline,
                                                                                     xDataTemporary[1]);
        printf("\n%g \t %d\n", temporaryInterpolantDerivative,i);

    }
    /*
    //Fill subarrays with data
    xDataTemporary[0] = xData[0];
    xDataTemporary[1] = xData[1];
    xDataTemporary[2] = xData[2];
    yDataTemporary[0] = yData[0];
    yDataTemporary[1] = yData[1];
    yDataTemporary[2] = yData[2];

    printf("\n%g\n", xDataTemporary[0]);
    printf("%g\n", xDataTemporary[1]);
    printf("%g\n", xDataTemporary[2]);

    printf("\n%g\n", yDataTemporary[0]);
    printf("%g\n", yDataTemporary[1]);
    printf("%g\n", yDataTemporary[2]);


    //First without for-loop to test my method
    quadSpline *temporaryQuadraticSpline = initialize_quadratic_spline(temporaryNumberOfPoints, xDataTemporary,
                                                                       yDataTemporary);

    double temporaryInterpolantDerivative = evaluate_quadratic_spline_derivative(temporaryQuadraticSpline,
                                                                                 xDataTemporary[1]);

    printf("\n%g\n", temporaryInterpolantDerivative);

    xDataTemporary[0] = xData[1];
    xDataTemporary[1] = xData[2];
    xDataTemporary[2] = xData[3];
    yDataTemporary[0] = yData[1];
    yDataTemporary[1] = yData[2];
    yDataTemporary[2] = yData[3];

    printf("\n%g\n", xDataTemporary[0]);
    printf("%g\n", xDataTemporary[1]);
    printf("%g\n", xDataTemporary[2]);

    printf("\n%g\n", yDataTemporary[0]);
    printf("%g\n", yDataTemporary[1]);
    printf("%g\n", yDataTemporary[2]);

    temporaryQuadraticSpline = initialize_quadratic_spline(temporaryNumberOfPoints, xDataTemporary, yDataTemporary);

    double temporaryInterpolantDerivative2 = evaluate_quadratic_spline_derivative(temporaryQuadraticSpline,
                                                                                  xDataTemporary[1]);

    printf("\n%g\n", temporaryInterpolantDerivative2);

     */


    //Settings for function
    double lowerLimit = xData[0];
    double upperLimit = 11;
    double absoluteError = 1e-6;
    double relativeError = 1e-6;
    size_t iterationLimit = 999;

    //Compute integral using GSL
    cosine_integral(numberOfPoints, xData, yData, lowerLimit, upperLimit, absoluteError, relativeError, iterationLimit);

    //Settings for GSL interpolationOld to compare
    gsl_interp *linearInterpolationGSL = gsl_interp_alloc(gsl_interp_linear, numberOfPoints);
    gsl_interp *quadraticInterpolationGSL = gsl_interp_alloc(gsl_interp_polynomial, numberOfPoints);
    gsl_interp *cubicInterpolationGSL = gsl_interp_alloc(gsl_interp_cspline, numberOfPoints);
    gsl_interp_init(linearInterpolationGSL, xData, yData, numberOfPoints);
    gsl_interp_init(quadraticInterpolationGSL, xData, yData, numberOfPoints);
    gsl_interp_init(cubicInterpolationGSL, xData, yData, numberOfPoints);


    //PART A
    printf("\nA: Linear spline interpolation \n\n");
    double resolution = fabs(xData[numberOfPoints - 1] - xData[0]) / numberOfSamples;

    for (double i = xData[0]; i < xData[numberOfPoints]; i += resolution)
    {

        double temporaryInterpolant = linear_spline(numberOfPoints, xData, yData, i);
        double temporaryInterpolantGSL = gsl_interp_eval(linearInterpolationGSL, xData, yData, i, NULL);
        double temporaryInterpolantIntegral = linear_spline_integrate(numberOfPoints, xData, yData, i);
        double temporaryInterpolantIntegralGSL = gsl_interp_eval_integ(linearInterpolationGSL, xData, yData, xData[0],
                                                                       i, NULL);

        //Output to linear data file
        fprintf(linearOutput, "%g\t%g\t%g\t%g\t%g\n", i, temporaryInterpolant, temporaryInterpolantGSL,
                temporaryInterpolantIntegral,
                temporaryInterpolantIntegralGSL);
    }
    printf("Results can be seen in linearSplinePlot.png\n\n");


    //PART B
    printf("B: Quadratic spline interpolation \n\n");

    quadSpline *quadraticSpline = initialize_quadratic_spline(numberOfPoints, xData, yData);

    for (double i = xData[0]; i < xData[numberOfPoints]; i += resolution)
    {
        double temporaryInterpolant = evaluate_quadratic_spline(quadraticSpline, i);
        double temporaryInterpolantGSL = gsl_interp_eval(quadraticInterpolationGSL, xData, yData, i, NULL);
        double temporaryInterpolantIntegral = evaluate_quadratic_spline_integral(quadraticSpline, i);
        double temporaryInterpolantIntegralGSL = gsl_interp_eval_integ(quadraticInterpolationGSL, xData, yData,
                                                                       xData[0], i, NULL);
        double temporaryInterpolantDerivative = evaluate_quadratic_spline_derivative(quadraticSpline, i);
        double temporaryInterpolantDerivativeGSL = gsl_interp_eval_deriv(quadraticInterpolationGSL, xData, yData, i,
                                                                         NULL);

        //Output to quadratic data file
        fprintf(quadraticOutput, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", i, temporaryInterpolant, temporaryInterpolantGSL,
                temporaryInterpolantIntegral, temporaryInterpolantIntegralGSL, temporaryInterpolantDerivative,
                temporaryInterpolantDerivativeGSL);
    }
    printf("Results can be seen in quadraticSplinePlot.png\n\n");

    //PART C
    printf("C: Cubic spline interpolation \n\n");

    cubicSpline *cubicSpline = initialize_cubic_spline(numberOfPoints, xData, yData);

    for (double i = xData[0]; i < xData[numberOfPoints]; i += resolution)
    {

        double temporaryInterpolant = evaluate_cubic_spline(cubicSpline, i);
        double temporaryInterpolantGSL = gsl_interp_eval(cubicInterpolationGSL, xData, yData, i, NULL);
        double temporaryInterpolantIntegral = evaluate_cubic_spline_integral(cubicSpline, i);
        double temporaryInterpolantIntegralGSL = gsl_interp_eval_integ(cubicInterpolationGSL, xData, yData, xData[0], i,
                                                                       NULL);
        double temporaryInterpolantDerivative = evaluate_cubic_spline_derivative(cubicSpline, i);
        double temporaryInterpolantDerivativeGSL = gsl_interp_eval_deriv(cubicInterpolationGSL, xData, yData, i, NULL);

        fprintf(cubicOutput, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", i, temporaryInterpolant, temporaryInterpolantGSL,
                temporaryInterpolantIntegral, temporaryInterpolantIntegralGSL, temporaryInterpolantDerivative,
                temporaryInterpolantDerivativeGSL);
    }
    printf("Results can be seen in cubicSplinePlot.png\n\n");

    //Free memory
    fclose(linearOutput);
    fclose(quadraticOutput);

    free_quadratic_spline(quadraticSpline);
    free_cubic_spline(cubicSpline);
    gsl_interp_free(linearInterpolationGSL);
    gsl_interp_free(quadraticInterpolationGSL);

    return 0;
}