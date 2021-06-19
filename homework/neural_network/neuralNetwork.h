#ifndef HAVE_NEURALNETWORK_H
#define HAVE_NEURALNETWORK_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct
{
    int numberOfNeurons;

    double (*target_function)(double);

    gsl_vector *paramters;
} neuralNetwork;

double random_number(unsigned int *seed);

double neural_network_response(neuralnetwork *neuralNetwork, double evaluationPoint);

void neural_network_train(neuralnetwork *neuralNetwork, gsl_vector *inputData, gsl_vector *labels);

neuralnetwork *neural_network_allocation(int numberOfNeuronsInput, double (*targetFunctionInput)(double));

void free_neural_network(neuralnetwork *neuralNetwork);


#endif //HAVE_NEURALNETWORK_H


