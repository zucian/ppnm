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

    //Insert data form file "cosData.txt" into arrays
    inputToArray(xData, yData, inputFile);

    //Esimate p_i of data
    estimate_derivative(numberOfPoints, xData, yData, pData);

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

    //Initializing sub spline for cosine function
    subSpline *subSplineCos = initialize_sub_spline(numberOfPoints, xData, yData, pData);

    //Compute sub spline value
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

    int numberOfSamplesJump = 500;
    double resolutionJump = fabs(xDataJump[numberOfPointsJump - 1] - xDataJump[0]) / numberOfSamplesJump;

    for (double i = xDataJump[0]; i < xDataJump[numberOfPointsJump-1]; i += resolutionJump)
    {
        double temporaryInterpolantSubSpline = evaluate_sub_spline(subSplineJump, i);
        fprintf(jumpOutput, "%g\t%g\n", i, temporaryInterpolantSubSpline);
    }

    //Free memory
    fclose(cubicOutput);
    fclose(subOutput);
    fclose(jumpOutput);

    free_cubic_spline(cubicSpline);
    free_sub_spline(subSplineCos);
    free_sub_spline(subSplineJump);

    return 0;
}