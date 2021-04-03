
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include "linearSpline.h"
#include "utilities.h"

double linear_spline(int nPoints, double* points, double* funcPoints, double evalPoints){

    int interval = binary_search(nPoints, points, evalPoints);
    double funcDifference = (funcPoints[interval+1] - funcPoints[interval]);
    double pointsDifference = (points[interval+1] - points[interval]);
    double slope = funcDifference/pointsDifference;

    double interpolationValue = funcPoints[interval] + slope*(evalPoints - points[interval]);
    return interpolationValue;
}

double linear_spline_integration(int nPoints, double* points, double* funcPoints, double evalPoints){

    int interval = binary_search(nPoints, points, evalPoints);

    double integral = 0;
    double funcDifference;
    double pointsDifference;
    double slope;

    for(int i=0; i <= interval; i++){
        funcDifference = (funcPoints[i+1] - funcPoints[i]);
        pointsDifference = (points[i+1] - points[i]);
        slope = funcDifference/pointsDifference;

        if(i >= interval){
            pointsDifference = (evalPoints - points[i]);
        }
        integral += funcPoints[i]*pointsDifference + slope*pointsDifference*pointsDifference/2;
    }
    return integral;
}