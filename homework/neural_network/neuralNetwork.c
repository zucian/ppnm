#include <time.h>
#include "neuralNetwork.h"
#include "minimization.h"
#include "utilities.h"


neuralNetwork *neural_network_allocation(int numberOfNeuronsInput, double (*targetFunctionInput)(double),
                                         double (*targetDerivativeInput)(double), double(*targetIntegralInput)(double))
{
    int numberOfParameters = 3;
    neuralNetwork *nn = (neuralNetwork *) malloc(sizeof(neuralNetwork));
    nn->parameters = gsl_vector_alloc(numberOfNeuronsInput * numberOfParameters);
    nn->targetFunction = targetFunctionInput;
    nn->targetDerivative = targetDerivativeInput;
    nn->targetIntegral = targetIntegralInput;
    nn->numberOfNeurons = numberOfNeuronsInput;

    return nn;
}

double neural_network_response(neuralNetwork *network, double evaluationPoint)
{
    int numberOfParameters = 2;
    int numberOfNeurons = network->numberOfNeurons;
    double response = 0;

    for (int i = 0; i < numberOfNeurons; ++i)
    {
        double neuronShift = gsl_vector_get(network->parameters, numberOfParameters * i);
        double neuronScale = gsl_vector_get(network->parameters, numberOfParameters * i + 1);
        double neuronEdge = gsl_vector_get(network->parameters, numberOfParameters * i + 2);

        response += (network->targetFunction((evaluationPoint - neuronShift) / neuronScale)) * neuronEdge;
    }
    return response;
}

double neural_network_response_derivative(neuralNetwork *network, double evaluationPoint)
{
    int numberOfParameters = 2;
    int numberOfNeurons = network->numberOfNeurons;
    double response = 0;

    for (int i = 0; i < numberOfNeurons; ++i)
    {
        double neuronShift = gsl_vector_get(network->parameters, numberOfParameters * i);
        double neuronScale = gsl_vector_get(network->parameters, numberOfParameters * i + 1);
        double neuronEdge = gsl_vector_get(network->parameters, numberOfParameters * i + 2);

        response += (network->targetDerivative((evaluationPoint - neuronShift) / neuronScale)) * neuronEdge /
                    neuronScale;
    }
    return response;
}

double neural_network_response_integral(neuralNetwork *network, double rightPoint, double leftPoint)
{
    int numberOfParameters = 2;
    int numberOfNeurons = network->numberOfNeurons;
    double response = 0;

    for (int i = 0; i < numberOfNeurons; ++i)
    {
        double neuronShift = gsl_vector_get(network->parameters, numberOfParameters * i);
        double neuronScale = gsl_vector_get(network->parameters, numberOfParameters * i + 1);
        double neuronEdge = gsl_vector_get(network->parameters, numberOfParameters * i + 2);

        response += (((network->targetIntegral)((leftPoint - neuronShift) / neuronScale)) * neuronEdge *
                     neuronScale) -
                    (((network->targetIntegral)((rightPoint - neuronShift) / neuronScale)) * neuronEdge *
                     neuronScale);
    }
    return response;
}

void neural_network_train(neuralNetwork *network, gsl_vector *inputData, gsl_vector *labels)
{
    unsigned int seed = time(NULL);
    int numberOfParameters = 3;
    int numberOfPoints = inputData->size;
    int numberOfNeurons = network->numberOfNeurons;

    double cost_function(gsl_vector *nextParameters)
    {
        int numberOfNeurons = network->numberOfNeurons;
        gsl_vector *updatedParameters = gsl_vector_alloc(numberOfParameters * numberOfNeurons);
        for (int i = 0; i < numberOfNeurons; ++i)
        {
            for (int j = 0; j < numberOfParameters; ++j)
            {
                gsl_vector_set(updatedParameters, numberOfParameters * i + j,
                               gsl_vector_get(nextParameters, numberOfParameters * i + j));
            }
        }

        network->parameters = updatedParameters; //Update parameters
        double cost = 0;

        for (int i = 0; i < numberOfPoints; ++i)
        {
            double evaluationPoint = gsl_vector_get(inputData, i);
            double pointLabel = gsl_vector_get(labels, i);
            double response = neural_network_response(network, evaluationPoint);
            cost += (response - pointLabel) * (response - pointLabel);
        }

        return cost;
    }

    double tolerance = 1e-5;
    gsl_vector *learnedParameters = gsl_vector_alloc(numberOfNeurons * numberOfParameters);

    double a = -5;
    double b = 5;

    for (int i = 0; i < numberOfNeurons; i++)
    {
        double weight = 1.000001;
        double scale = 1;
        double shift = -5+(b-a)*i/(numberOfNeurons-1);

        gsl_vector_set(learnedParameters,3*i,1);
        gsl_vector_set(learnedParameters, 3*i+1,scale);
        gsl_vector_set(learnedParameters,3*i+2,weight);
    }

    quasi_newton_method(cost_function, learnedParameters, tolerance);

    network->parameters = learnedParameters;

}


void free_neural_network(neuralNetwork *network)
{
    gsl_vector_free(network->parameters);
    free(network);
}


