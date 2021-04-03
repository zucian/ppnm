#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include "utilities.h"
#include "linearSpline.h"

int main(){

    int numberOfPoints = 20;
    int numberOfSamples = 500;

    FILE* linearSpline = fopen("out.linearSpline.txt","w");
    FILE* sineData = fopen("out.sineData.txt","w");

    for(int i=0;i<numberOfPoints;i++){
        fprintf(sineData,"%g\t%g\n",(double) i/2., sin((double) i/2.));
    }

    double* xData = malloc(numberOfPoints*sizeof(double));
    double* yData = malloc(numberOfPoints*sizeof(double));

    input_to_array(xData, yData, "sineData");

    double lowerLimit = xData[0];
    double upperLimit = 9.5;
    double absoluteError = 1e-6;
    double relativeError = 1e-6;
    size_t iterationLimit = 1000;

    //Exercise A
    gsl_interp* gslInterpolationLinear = gsl_interp_alloc(gsl_interp_linear, numberOfPoints);
    gsl_interp_init(gslInterpolationLinear, xData, yData, numberOfPoints);

    double resolution = fabs(xData[numberOfPoints-1]-xData[0])/numberOfSamples;

    for(double evalPoint = xData[0]; evalPoint < xData[numberOfPoints]; evalPoint += resolution){

        double interpolationPart = linear_spline(numberOfPoints,xData,yData,evalPoint);
        double interpolationPartGsl = gsl_interp_eval(gslInterpolationLinear, xData, yData, evalPoint, NULL);
        double interpolationIntegralPart = linear_spline_integration(numberOfPoints,xData,yData,evalPoint);
        double interpolationIntegralPartGsl = gsl_interp_eval_integ(gslInterpolationLinear,xData,Ydata,xData[0],evalPoint, NULL);

        fprintf(linearSpline,"%g\t%g\t%g\t%g\t%g\n", evalPoint,interpolationPart,interpolationPartGsl,interpolationIntegralPart,interpolationIntegralPartGsl);
    }

    return 0;
}