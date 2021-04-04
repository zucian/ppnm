#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include "utilities.h"
#include "linearSpline.h"
#include "quadraticSpline.h"

int main() {

    int numberOfPoints = 20;
    int numberOfSamples = 500;

    FILE *sineData = fopen("out.sineData.txt", "w");

    for (int i = 0; i < numberOfPoints; i++) {
        fprintf(sineData, "%g\t%g\n", (double) i / 2., sin((double) i / 2.));
    }

    double *xData = malloc(numberOfPoints * sizeof(double));
    double *yData = malloc(numberOfPoints * sizeof(double));

    input_to_array(xData, yData, "sineData");

    double lowerLimit = xData[0];
    double upperLimit = 9.5;
    double absoluteError = 5e-5;
    double relativeError = 5e-5;
    size_t iterationLimit = 1000;

    //Exercise A
    FILE *linearSpline = fopen("out.linearSpline.txt", "w");
    gsl_interp *gslInterpolationLinear = gsl_interp_alloc(gsl_interp_linear, numberOfPoints);
    gsl_interp_init(gslInterpolationLinear, xData, yData, numberOfPoints);

    double resolution = fabs(xData[numberOfPoints - 1] - xData[0]) / numberOfSamples;

    for (double pointsToEvaluateInterpolant = xData[0]; pointsToEvaluateInterpolant < xData[numberOfPoints]; pointsToEvaluateInterpolant += resolution)
    {
        double interpolationPart = linear_spline(numberOfPoints, xData, yData, pointsToEvaluateInterpolant);
        double interpolationIntegralPart = linear_spline_integration(numberOfPoints, xData, yData, pointsToEvaluateInterpolant);
        double interpolationPartGsl = gsl_interp_eval(gslInterpolationLinear, xData, yData, pointsToEvaluateInterpolant, NULL);
        double interpolationIntegralPartGsl = gsl_interp_eval_integ(gslInterpolationLinear, xData, Ydata, xData[0],
                                                                    pointsToEvaluateInterpolant, NULL);

        fprintf(linearSpline, "%g\t%g\t%g\t%g\t%g\n", pointsToEvaluateInterpolant, interpolationPart, interpolationIntegralPart, interpolationPartGsl,
                 interpolationIntegralPartGsl);
    }

    //Exercise B
    FILE *quadraticSpline = fopen("out.quadraticSpline.txt", "w");
    gsl_interp *gslInterpolationQuadratic = gsl_interp_alloc(gsl_interp_polynomial, numberOfPoints);
    gsl_interp_init(gslInterpolationQuadratic, xData, yData, numberOfPoints);

    quadSpline *quadraticSpline = initialize_quadratic_spline(numberOfPoints, xData, yData);

    for(double pointsToEvaluateInterpolant = xData[0]; pointsToEvaluateInterpolant < xData[numberOfPoints]; pointsToEvaluateInterpolant += resolution)
    {
        double interpolationPart = evaluate_quadratic_spline(quadraticSpline, pointsToEvaluateInterpolant);
        double interpolationIntegralPart = evaluate_quadratic_spline_integral(quadraticSpline,pointsToEvaluateInterpolant);
        double interpolationDerivativePart = evaluate_quadratic_spline_derivative(quadraticSpline,pointsToEvaluateInterpolant);
        double interpolationPartGsl = gsl_interp_eval(gslInterpolationQuadratic, xData, yData, pointsToEvaluateInterpolant, NULL);
        double interpolationIntegralPartGsl = gsl_interp_eval_integ(gslInterpolationQuadratic, xData, yData, xData[0], pointsToEvaluateInterpolant, NULL);
        double interpolationDerivativePartGsl = gsl_interp_eval_deriv(gslInterpolationQuadratic, xData, yData, pointsToEvaluateInterpolant, NULL);

        fprintf(linearSpline, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", pointsToEvaluateInterpolant, interpolationPart, interpolationIntegralPart, interpolationDerivativePart interpolationPartGsl,
                interpolationIntegralPartGsl, interpolationDerivativePartGsl);
    }

    return 0;
}