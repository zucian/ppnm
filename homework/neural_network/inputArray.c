#include <stdio.h>
#include <stdlib.h>
#include "inputArray.h"

void inputArray(double *xData, double *yData, double *zData, char *inputFilename)
{
    int numberOfReturnValues = 3;
    int input = numberOfReturnValues;
    double argumentX;
    double argumentY;
    double argumentZ;

    FILE *inputFileStream = fopen(inputFilename, "r");

    int i = 0;
    while (input != EOF)//end of file = EOF
    {
        if (input == numberOfReturnValues)
        {
            input = fscanf(inputFileStream, "%lg \t %lg \t %lg", &argumentX, &argumentY, &argumentZ);
            xData[i] = argumentX;
            yData[i] = argumentY;
            zData[i] = argumentZ;

            i++;
        }
        else
        {
            fprintf(stderr,"Failed to read input \n");
            exit(-1);
        }
    }
    fclose(inputFileStream);
}