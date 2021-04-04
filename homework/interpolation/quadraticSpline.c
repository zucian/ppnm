#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "utilities.h"
#include "quadraticSpline.h"

quadSpline* initialize_quadratic_spline(int numberOfPoints, double* points, double* functionValueOfPoints){

    quadSpline* spline = (quadSpline*) malloc(sizeof(quadSpline));

    int numberOfEquations = numberOfPoints-1;

    spline->numberOfPoins = numberOfPoints;
    spline->points = (double*) (numberOfEquations*sizeof(double));
    spline->functionValueOfPoints = (double*) (numberOfEquations*sizeof(double));
    spline->firstCoefficient = (double*) malloc(numberOfEquations*sizeof(double));
    spline->secondCoefficient = (double*) malloc(numberOfEquations*sizeof(double));

    for(int i=0; i<numberOfPoints; i++){
        spline->points[i] = points[i];
        spline->functionValueOfPoints[i] = functionValueOfPoints[i];
    }

    double pointsDifference[numberOfPoints-1];
    double slope[numberOfPoints-1];

    for(int i=0; i<numberOfPoints-1; i++){
        pointsDifference[i] = points[i+1]-points[i];
        slope[i] = (functionValueOfPoints[i+1]-functionValueOfPoints[i])/pointsDifference[i];
    }

    //Forward recursion to compute second coefficient
    spline->secondCoefficient[0] = 0;
    for(int i=0; i<numberOfPoints-2; i++){
        spline->secondCoefficient[i+1] = (slope[i+1]-slope[i]-(spline->secondCoefficient[i])*pointsDifference[i])/pointsDifference[i+1];
    }

    //Backward recursion to compute second coefficient
    spline->secondCoefficient[numberOfPoints-2] /= 2;
    for(int i=numberOfPoints-3; i>=0; i--){
        spline->secondCoefficient[i] = (slope[i+1]-slope[i]-(spline->secondCoefficient[i+1])*pointsDifference[i+1])/pointsDifference[i];
    }

    //First coefficient
    for(int i=0; i<numberOfPoints-1; i++){
        spline->firstCoefficient[i] = slope[i]-(spline->secondCoefficient[i])*pointsDifference[i];
    }

    return spline;
}

double evaluate_quadratic_spline(quadSpline* spline, double pointsToEvaluateInterpolant){

    int whichInterval = binary_search(spline->numberOfPoints, spline->points, pointsToEvaluateInterpolant);
    double pointsDifference = pointsToEvaluateInterpolant-(spline->points[whichInterval]);
    double interpolantValue = (spline->functionValueOfPoints[whichInterval])+pointsDifference*((spline->firstCoefficient[whichInterval])+pointsDifference*(secondCoefficient[whichInterval]));
    return interpolantValue;
}

double evaluate_quadratic_spline_derivative(quadSpline* spline, double pointsToEvaluateInterpolant){

    int whichInterval = binary_search(spline->numberOfPoints, spline->points, pointsToEvaluateInterpolant);
    double pointsDifference = pointsToEvaluateInterpolant-(spline->points[whichInterval]);
}