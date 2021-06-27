#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include "utilities.h"
#include "integrateFunction.h"
#include "quadraticSpline.h"
#include "cubicSpline.h"
#include "subSpline.h"

double gsl_cos(double x, void *params)
{
    return cos(x);
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
    FILE *cubicOutput = fopen(argv[2], "w");
    FILE *subOutput = fopen(argv[3], "w");
    FILE *jumpOutput = fopen(argv[4], "w");
    FILE *jumpData = fopen(argv[5], "w");

    double *xData = malloc(numberOfPoints * sizeof(double));
    double *yData = malloc(numberOfPoints * sizeof(double));
    double *pData = malloc(numberOfPoints * sizeof(double));

    inputToArray(xData, yData, inputFile);

    estimate_derivative(numberOfPoints, xData, yData, pData);

    /*
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
    int counter = 1;
    for (int i = 1; i < numberOfPoints - 1; i++)
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

        //If this is first point, use polynomial to estimate p_1
        if (counter == 1)
        {
            double temporaryInterpolantDerivativeStart = evaluate_quadratic_spline_derivative(temporaryQuadraticSpline,
                                                                                              xDataTemporary[0]);
            //printf("\n%g \t %d\n", temporaryInterpolantDerivativeStart, i);
            //fprintf(dataFile,"%g\t%g\t%g\n", xData[i - 1], yData[i - 1],temporaryInterpolantDerivativeStart);
            pData[i - 1] = temporaryInterpolantDerivativeStart;
        }

        //Evaluate p_i
        double temporaryInterpolantDerivative = evaluate_quadratic_spline_derivative(temporaryQuadraticSpline,
                                                                                     xDataTemporary[1]);
        pData[i] = temporaryInterpolantDerivative;
        //printf("\n%g \t %d\n", temporaryInterpolantDerivative, i + 1);
        //fprintf(dataFile,"%g\t%g\t%g\n", xData[i], yData[i],temporaryInterpolantDerivative);

        //Add one to counter that checks if first or last point
        counter++;

        //If this is last point, use polynomial to estimate p_n
        if (counter == numberOfPoints - 1)
        {
            double temporaryInterpolantDerivativeEnd = evaluate_quadratic_spline_derivative(temporaryQuadraticSpline,
                                                                                            xDataTemporary[2]);
            //printf("\n%g \t %d\n", temporaryInterpolantDerivativeEnd, i + 2);
            //fprintf(dataFile,"%g\t%g\t%g\n", xData[i + 1], yData[i + 1],temporaryInterpolantDerivativeEnd);
            pData[i + 1] = temporaryInterpolantDerivativeEnd;
        }

    }
    */

    //Spacing between data points
    double resolution = fabs(xData[numberOfPoints - 1] - xData[0]) / numberOfSamples;

    //Doing cubic spline again to have something to compare to
    cubicSpline *cubicSpline = initialize_cubic_spline(numberOfPoints, xData, yData);

    for (double i = xData[0]; i < xData[numberOfPoints]; i += resolution)
    {
        double temporaryInterpolant = evaluate_cubic_spline(cubicSpline, i);
        double temporaryInterpolantIntegral = evaluate_cubic_spline_integral(cubicSpline, i);
        double temporaryInterpolantDerivative = evaluate_cubic_spline_derivative(cubicSpline, i);

        fprintf(cubicOutput, "%g\t%g\t%g\t%g\n", i, temporaryInterpolant,
                temporaryInterpolantIntegral, temporaryInterpolantDerivative
        );
    }


    //EXAM PART

    //Print loop to show xData, yData, and pData content
    /*
    for (int i = 0; i < numberOfPoints; i++)
    {
        printf("%d\t%g\t%g\t%g\n", i, xData[i], yData[i], pData[i]);
    }
     */


    subSpline *subSplineCos = initialize_sub_spline(numberOfPoints, xData, yData, pData);

    for (double i = xData[0]; i < xData[numberOfPoints]; i += resolution)
    {
        double temporaryInterpolantSubSpline = evaluate_sub_spline(subSplineCos, i);
        fprintf(subOutput, "%g\t%g\n", i, temporaryInterpolantSubSpline);
    }

    //Test on some new homemade data with sudden jump
    int numberOfPointsJump = 8;
    double xDataJump[8] = {0.0, 0.5, 1.0, 1.55, 1.95, 2.5, 3.0, 3.5};
    double yDataJump[8] = {0.23, 0.35, 0.25, 0.33, 1.37, 1.44, 1.25, 1.51};
    double *pDataJump = malloc(numberOfPointsJump * sizeof(double));

    int counter = 0;
    for (int i = 0; i < numberOfPointsJump; i++)
    {
        counter++;
        if (counter == numberOfPointsJump)
        {
            fprintf(jumpData,"%g\t%g", xDataJump[i],yDataJump[i]);
        }
        else
        {
            fprintf(jumpData,"%g\t%g\n", xDataJump[i],yDataJump[i]);
        }
    }

    estimate_derivative(numberOfPointsJump,xDataJump,yDataJump,pDataJump);

    subSpline *subSplineJump = initialize_sub_spline(numberOfPointsJump,xDataJump,yDataJump,pDataJump);

    double resolutionJump = fabs(xDataJump[numberOfPointsJump - 1] - xDataJump[0]) / numberOfSamples;

    printf("%g\n", xDataJump[0]);
    printf("%g\n", xDataJump[numberOfPointsJump-1]);
    for (double i = xDataJump[0]; i < xDataJump[numberOfPointsJump-1]; i += resolutionJump)
    {
        double temporaryInterpolantSubSpline = evaluate_sub_spline(subSplineJump, i);
        fprintf(jumpOutput, "%g\t%g\n", i, temporaryInterpolantSubSpline);
    }


    //Free memory
    fclose(cubicOutput);
    fclose(subOutput);

    free_cubic_spline(cubicSpline);
    free_sub_spline(subSplineCos);

    return 0;
}