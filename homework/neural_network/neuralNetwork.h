#ifndef HAVE_NEURALNETWORK_H
#define HAVE_NEURALNETWORK_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct
{
    int numberOfNeurons;

    double (*targetFunction)(double);

    double (*targetDerivative)(double);

    double (*targetIntegral)(double);

    gsl_vector *parameters;
} neuralNetwork;

double neural_network_response(neuralNetwork *neuralNetwork, double evaluationPoint);

void neural_network_train(neuralNetwork *neuralNetwork, gsl_vector *inputData, gsl_vector *labels);

double neural_network_response_derivative(neuralNetwork *neuralNetwork, double evaluationPoint);

double neural_network_response_integral(neuralNetwork *network, double rightPoint, double leftPoint);

neuralNetwork *neural_network_allocation(int numberOfNeuronsInput, double (*targetFunctionInput)(double),
                                         double (*targetDerivativeInput)(double), double(*targetIntegralInput)(double));

void free_neural_network(neuralNetwork *neuralNetwork);


#endif //HAVE_NEURALNETWORK_H


