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
        fprintf(sineData,"%g\t%g\n",(double) i/2., sin((double) i/2));
    }


    double* xData = malloc(numberOfPoints*sizeof(double));
    double* yData = malloc(numberOfPoints*sizeof(double));

    input_to_array(xData, yData, "sineData");

    printf(xData);
    printf(yData);
    return 0;
}