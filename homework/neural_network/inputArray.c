#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "inputArray.h"

void input_to_array(int numberOfDataPoints, gsl_vector *xData, gsl_vector *yData, char *inputFilename)
{
    int numberOfReturnValues = 2;
    int input = numberOfReturnValues;
    double argumentX;
    double argumentY;

    FILE *inputFileStream = fopen(inputFilename, "r");

    int i = 0;
    while (input != EOF && i < numberOfDataPoints)//end of file = EOF
    {
        if (input == numberOfReturnValues)
        {
            input = fscanf(inputFileStream, "%lg \t %lg", &argumentX, &argumentY);
            gsl_vector_set(xData, i, argumentX);
            gsl_vector_set(yData, i, argumentY);

            i++;
        }
        else
        {
            fprintf(stderr, "Failed to read input \n");
            exit(-1);
        }
    }
    fclose(inputFileStream);
}