#include "differentialEquations.h"
#include <math.h>
#include <gsl/gsl_vector.h>

void harmonic_function(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative)
{
    gsl_vector_set(functionDerivative, 0, gsl_vector_get(functionValue, 1));
    gsl_vector_set(functionDerivative, 1, -gsl_vector_get(functionValue, 0));
}

void SIR_model(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative)
{
    double population = 5806000; //population of denmark 2019
    double recoveryTime = 21; //Three weeks recovery
    double contactTime = 3; //Time between contacts
    double susceptible = gsl_vector_get(functionValue, 0);
    double infected = gsl_vector_get(functionValue, 1);
    double derivativeSusceptible = -(susceptible * infected) / (contactTime * population);
    double derivativeRemoved = infected / recoveryTime;
    double derivativeInfected = -derivativeSusceptible - derivativeRemoved;

    gsl_vector_set(functionDerivative, 0, derivativeSusceptible);
    gsl_vector_set(functionDerivative, 1, derivativeInfected);
    gsl_vector_set(functionDerivative, 2, derivativeRemoved);
}

void SIR_model_new_contact_time(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative)
{
    double population = 5806000; //population of denmark 2019
    double recoveryTime = 21; //Three weeks recovery
    double contactTime = 6; //Time between contacts, doubled for this equation
    double susceptible = gsl_vector_get(functionValue, 0);
    double infected = gsl_vector_get(functionValue, 1);
    double derivativeSusceptible = -(susceptible * infected) / (contactTime * population);
    double derivativeRemoved = infected / recoveryTime;
    double derivativeInfected = -derivativeSusceptible - derivativeRemoved;

    gsl_vector_set(functionDerivative, 0, derivativeSusceptible);
    gsl_vector_set(functionDerivative, 1, derivativeInfected);
    gsl_vector_set(functionDerivative, 2, derivativeRemoved);
}

void three_body_problem(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative)
{
    double mass1 = 1;
    double mass2 = 1;
    double mass3 = 1;
    double x1 = gsl_vector_get(functionValue, 0);
    double x2 = gsl_vector_get(functionValue, 2);
    double x3 = gsl_vector_get(functionValue, 4);
    double y1 = gsl_vector_get(functionValue, 1);
    double y2 = gsl_vector_get(functionValue, 3);
    double y3 = gsl_vector_get(functionValue, 5);
    double G = 1;
    gsl_vector_set(functionDerivative, 0, gsl_vector_get(functionValue, 6));
    gsl_vector_set(functionDerivative, 1, gsl_vector_get(functionValue, 7));
    gsl_vector_set(functionDerivative, 2, gsl_vector_get(functionValue, 8));
    gsl_vector_set(functionDerivative, 3, gsl_vector_get(functionValue, 9));
    gsl_vector_set(functionDerivative, 4, gsl_vector_get(functionValue, 10));
    gsl_vector_set(functionDerivative, 5, gsl_vector_get(functionValue, 11));
    double distance12 = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
    double distance13 = sqrt(pow((x3 - x1), 2) + pow((y3 - y1), 2));
    double distance23 = sqrt(pow((x2 - x3), 2) + pow((y2 - y3), 2));
    gsl_vector_set(functionDerivative, 6, (x2 - x1) / pow(distance12, 3) + (x3 - x1) / pow(distance13, 3));
    gsl_vector_set(functionDerivative, 7, (y2 - y1) / pow(distance12, 3) + (y3 - y1) / pow(distance13, 3));
    gsl_vector_set(functionDerivative, 8, (x3 - x2) / pow(distance23, 3) + (x1 - x2) / pow(distance12, 3));
    gsl_vector_set(functionDerivative, 9, (y3 - y2) / pow(distance23, 3) + (y1 - y2) / pow(distance12, 3));
    gsl_vector_set(functionDerivative, 10, (x2 - x3) / pow(distance23, 3) + (x1 - x3) / pow(distance13, 3));
    gsl_vector_set(functionDerivative, 11, (y2 - y3) / pow(distance23, 3) + (y1 - y3) / pow(distance13, 3));
}