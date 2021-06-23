#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "neuralNetwork.h"
#include "inputArray.h"

double negative_sin(double x)
{
    return -sin(x); //Need function for derivative of cosine
}

int main(int argc, char *argv[])
{
    //PART A
    printf("A: Construct simple artificial neural network\n\n");

    int numberOfNeurons = 5;
    int numberOfDataPoints = 20;

    //Load cosine data
    gsl_vector *xData = gsl_vector_alloc(numberOfDataPoints);
    gsl_vector *yData = gsl_vector_alloc(numberOfDataPoints);
    input_to_array(numberOfDataPoints, xData, yData, argv[1]);

    printf("Initialize neural network with %d neurons, one hidden layer \n\n", numberOfNeurons);
    neuralNetwork *network = neural_network_allocation(numberOfNeurons, &cos, &negative_sin, &sin);

    printf("Training network: \n\n");
    neural_network_train(network, xData, yData);

    FILE *outputFile = fopen(argv[2], "w");

    double lowerLimit = 0;
    double upperLimit = 11;

    for (double i = lowerLimit; i <= upperLimit; i += 1.0 / 8)
    {
        fprintf(outputFile, "%10g %10g %10g %10g %10g %10g %10g\n", i, neural_network_response(network, i), cos(i),
                neural_network_response_derivative(network, i), -sin(i),
                neural_network_response_integral(network, 0, i), sin(i));
    }

    printf("Part A and B\n\n");
    printf("Network predicted form of cos(x) and its derivative + antiderivative can be seen in file networkPrediction.png\n\n");


    fclose(outputFile);

    free_neural_network(network);

    return 0;
}