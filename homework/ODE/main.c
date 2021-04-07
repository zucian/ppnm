#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void matrix_print(gsl_matrix *A, FILE *file)
{
    for (int i = 0; i < (A->size1); i++)
    {
        for (int j = 0; j < (A->size2); j++)
        {
            double Aij = gsl_matrix_get(A, i, j);
            fprintf(file, "%0.3g ", Aij);
        }
        fprintf(file, "\n");
    }
}

void rkstep12(void f(double x, gsl_vector *y, gsl_vector *derivativeOfY), double x, gsl_vector *yEvaluatedAtX,
              double stepSize, gsl_vector *yEvaluatedAtXWithStep, gsl_vector *errorEstimate)
{
    int dimensionOfYRange = yEvaluatedAtX->size;

    //Allocating memory for vectors
    gsl_vector *k0 = gsl_vector_alloc(dimensionOfYRange);
    gsl_vector *k1 = gsl_vector_alloc(dimensionOfYRange);
    gsl_vector *yFirstOrder = gsl_vector_alloc(dimensionOfYRange);

    //First order
    f(x, yEvaluatedAtX, k0);
    for (int i = 0; i < dimensionOfYRange; i++)
    {
        double k0Part = gsl_vector_get(k0, i);
        double yEvaluatedAtXPart = gsl_vector_get(yEvaluatedAtX, i);
        double yFirstOrderPart = yEvaluatedAtXPart + stepSize / 2 * k0part;

        gsl_vector_set(yFirstOrder, i, yFirstOrderPart);
    }

    //Second order
    f(x + stepSize / 2, yFirstOrder, k1);
    for (int i = 0; i < dimensionOfYRange; i++)
    {
        double k1part = gsl_vector_get(k1, i);
        double yEvaluatedAtXPart = gsl_vector_get(yEvaluatedAtX, i);
        double ySecondOrderPart = yEvaluatedAtXPart + h * k1part;

        gsl_vector_set(yEvaluatedAtXWithStep, i, ySecondOrderPart)
    }

    //Estimation of error
    for (int i = 0; i < dimensionOfYRange; i++)
    {
        double k0part = gsl_vector_get(k0, i);
        double k1part = gsl_vector_get(k1, i);
        double errorPart = stepSize / 2 * (k0part - k1part);

        gsl_vector_set(errorEstimate, i, errorPart);
    }

    //Freeing allocated memory
    gsl_vector_free(k0);
    gsl_vector_free(k1);
    gsl_vector_free(yFirstOrder);
}

void SIR_model_denmark(double x, gsl_vector *y, gsl_vector *derivativeOfY)
{
    /*
    vector y contains susceptible population in first entrance, infected in second, removed in thirds
    */

    //Parameters for model
    double populationOfDenmark = 5.8e6; //Population of denmark 2019
    double contactsOfInfected = 1.8;
    double recoveryTime = 14.; //Unit of days
    double timeBetweenContacts = recoveryTime / contactsOfInfected;

    //The three compartments of the population in the SIR model
    double susceptible = gsl_vector_get(y, 0);
    double infected = gsl_vector_get(y, 1);
    double removed = gsl_vector_get(y, 2);

    double derivativeOfSusceptible = -infected * susceptible / (populationOfDenmark * timeBetweenContacts);
    double derivativeOfInfected =
            infected * susceptible / (populationOfDenmark * timeBetweenContacts) - infected / recoveryTime;
    double derivativeOfRemoved = infected / recoveryTime

    gsl_vector_set(derivativeOfY, 0, derivativeOfSusceptible);
    gsl_vector_set(derivativeOfY, 1, derivativeOfInfected);
    gsl_vector_set(derivativeOfY, 2, derivativeOfRemoved);
}