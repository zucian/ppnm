#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "neuralNetwork.h"
#include "inputArray.h"

void vector_print(char *string, gsl_vector *vector)
{
    printf("%s\n", string);
    for (int iter = 0; iter < vector->size; iter++)
    {
        printf("%10g ", gsl_vector_get(vector, iter));
    }
    printf("\n");
}

void print_matrix(int numOfRows, gsl_matrix *matrixToPrint, char *string)
{
    printf("\n%s\n", string);
    for (int rowId = 0; rowId < numOfRows; rowId++)
    {
        gsl_vector_view matrixToPrint_row = gsl_matrix_row(matrixToPrint, rowId);
        gsl_vector *vector = &matrixToPrint_row.vector;
        for (int iter = 0; iter < vector->size; iter++)
        {
            if (gsl_vector_get(vector, iter) > 1e-10)
            {
                printf("%10g\t", gsl_vector_get(vector, iter));
            }
            else
            { printf("%10g\t", 0.0); }
        }
        printf("\n");
    }
}

void set_data_symmetric(gsl_matrix *inputMatrix, unsigned int *seed)
{
    for (int rowId = 0; rowId < (inputMatrix->size1); rowId++)
    {
        gsl_matrix_set(inputMatrix, rowId, rowId, randomNumber(seed));

        for (int colId = rowId + 1; colId < inputMatrix->size2; colId++)
        {
            double matrixElement = randomNumber(seed);
            gsl_matrix_set(inputMatrix, rowId, colId, matrixElement);
            gsl_matrix_set(inputMatrix, colId, rowId, matrixElement);
        }
    }
}

int main()
{
    //PART A
    printf("A: Construct simple artificial neural network\n\n");

    int numberOfNeurons = 5;
    int numberOfDataPoints = 20;

    //Load cosine data
    gsl_vector *xData = gsl_vector_alloc(numberOfDataPoints);
    gsl_vector *yData = gsl_vector_alloc(numberOfDataPoints);
    inputArray(numberOfDataPoints, xData, yData, "cosData.txt");

    printf("Initialize neural network with %d neurons, one hidden layer \n",numberOfNeurons);
    neuralNetwork *network = neural_network_allocation(numberOfNeurons, &sin);

    printf("Training network: \n");
    neural_network_train(network, xData, yData);

    FILE* outputFileStream = fopen()

    printf()



    return 0;
}