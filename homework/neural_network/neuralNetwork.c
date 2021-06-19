#include <time.h>
#include "neuralNetwork.h"
#include "minimization.h"
#include "utilities.h"

double random_number(unsigned int *seed)
{
    double maxRand = (double) RAND_MAX;           // Maximum random number, cast to double
    double randNum = (double) rand_r(seed);     // Generate pseudo-random number from seed, cast to double
    return randNum / maxRand;
}

double neural_network_response(neuralnetwork *neuralNetwork, double evaluationPoint)
{
    int numberOfParameters = 2;
    int numberOfNeurons = neuralNetwork->numberOfNeurons;
    double response = 0;

    for (int i = 0; i < numberOfNeurons; i++)
    {
        double neuronShift = gsl_vector_get(neuralNetwork->parameters, numberOfParameters * i);
        double neuronScale = gsl_vector_get(neuralNetwork->parameters, numberOfParameters * i + 1);
        double neuronEdge = gsl_vector_get(neuralNetwork->parameters, numberOfParameters * i + 2);

        response += (neuralNetwork->target_function((evaluationPoint - neuronShift) / neuronScale)) * neuronEdge;
    }

    return response;
}

void neural_network_train(neuralnetwork *neuralNetwork, gsl_vector *inputData, gsl_vector *labels)
{
    unsigned int seed = time(NULL);
    int numberOfParameters = 3;
    int numberOfPoints = inputData->size;
    int numberOfNeurons = neuralNetwork->numberOfNeurons;

    double cost_function(gsl_vector *nextParameters)
    {
        int numberOfNeurons = neuralNetwork->numberOfNeurons;
        gsl_vector *updatedParameters = gsl_vector_alloc(numberofParameters * numberOfNeurons);

        for (int i = 0; i < numberOfNeurons; i++)
        {
            for (int j = 0; j < numberofParameters; j++)
            {
                gsl_vector_set(updatedParameters, numberOfParameters * i + j,
                               gsl_vector_get(nextParameters, numberOfParameters * i + j));
            }
        }

        neuralNetwork->parameters = updatedParameters; //Update parameters
        double cost = 0;

        for (int i = 0; i < numberOfPoints; i++)
        {
            double evaluationPoint = gsl_vector_get(inputData, i);
            double pointLabel = gsl_vector_get(labels, i);
            double response = neural_network_response(neuralNetwork, evaluationPoint);
            cost += (response - pointLabel) * (response - pointLabel);
        }

        return cost;
    }

    double tolerance = 1e-5;
    gsl_vector *learnedParamters = gsl_vector_alloc(numberOfNeurons * numberOfParameters);
    for (int i = 0; i < numberOfNeurons * numberOfParameters; i++)
    {
        gsl_vector_set(learnedParameters, i, random_number(&seed));
    }
    quasi_newton_method(cost_function, learnedParameters, tolerance);

    neuralNetwork->paramters = learnedParameters;

}

neuralnetwork *neural_network_allocation(int numberOfNeuronsInput, double (*targetFunctionInput)(double))
{
    int numberOfParameters = 3;
    neuralnetwork *neuralNetwork = (neuralnetwork *) malloc(sizeof(neuralnetwork));
    neuralNetwork->parameters = gsl_vector_alloc(numberOfNeurons * numberOfParameters);
    neuralNetwork->targetFunction = targetFunctionInput;
    neuralNetwork->numberOfNeurons = numberOfNeuronsInput;

    return neuralNetwork;
}

void free_neural_network(neuralnetwork *neuralNetwork)
{
    gsl_vector_free(neuralNetwork->parameters);
    free(neuralNetwork);
}


