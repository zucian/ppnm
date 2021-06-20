#include <math.h>
#include "rungeKutta.h"
#include "differentialEquations.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void
rk_step12(void (*function)(double, gsl_vector *, gsl_vector *), double variable, gsl_vector *functionValue, double step,
          gsl_vector *functionStep, gsl_vector *error)
{
    int order = functionValue->size; //Order of differential equation
    gsl_vector *tangent0 = gsl_vector_alloc(order); //k_0
    gsl_vector *tangent12 = gsl_vector_alloc(order); //k_1/2
    gsl_vector *temporaryFunctionValue = gsl_vector_alloc(order);

    function(variable, functionValue, tangent0); //RHS of ODE

    for (int i = 0; i < order; i++)
    {
        double functionValueI = gsl_vector_get(functionValue, i);
        double tangent0I = gsl_vector_get(tangent0, i);
        double temporaryFunctionValueI = functionValueI + 0.5 * tangent0I * step;
        gsl_vector_set(temporaryFunctionValue, i, temporaryFunctionValueI);
    }

    function(variable + step * 0.5, temporaryFunctionValue, tangent12);

    //Advance solution
    for (int i = 0; i < order; i++)
    {
        double functionValueI = gsl_vector_get(functionValue, i);
        double tangent12I = gsl_vector_get(tangent12, i);
        double temporaryFunctionStepI = functionValueI + tangent12I * step;
        gsl_vector_set(functionStep, i, temporaryFunctionStepI);
    }

    //Compute error estimate, method from book
    for (int i = 0; i < order; i++)
    {
        double tangent0I = gsl_vector_get(tangent0, i);
        double tangent12I = gsl_vector_get(tangent12, i);
        double temporaryErrorI = (tangent0I - tangent12I) * step / 2;
        gsl_vector_set(error, i, temporaryErrorI);
    }

    //Free memory
    gsl_vector_free(tangent0);
    gsl_vector_free(tangent12);
    gsl_vector_free(temporaryFunctionValue);
}

void rk_driver(void (*function)(double, gsl_vector *, gsl_vector *), double leftEndpoint, gsl_vector *functionValueLeft,
               double rightEndPoint, gsl_vector *functionValueRight, double step, double absoluteAccuracy,
               double relativeAccuracy, FILE *pathToFile)
{
    int order = functionValueLeft->size; //Order of differential equation
    double tolerance;
    double error;
    double normFunctionValue;
    gsl_vector *functionStep = gsl_vector_alloc(order);
    gsl_vector *functionError = gsl_vector_alloc(order);
    gsl_vector *currentFunctionValue = gsl_vector_alloc(order);
    gsl_vector_memcpy(currentFunctionValue, functionValueLeft);
    double position = leftEndpoint;

    while (position < rightEndPoint)//Go through entire interval
    {
        if (pathToFile != NULL)
        {
            fprintf(pathToFile, "%.5g \t", position);

            for (int i = 0; i < order; i++)
            {
                fprintf(pathToFile, "%.5g \t", gsl_vector_get(currentFunctionValue, i));
            }
            if (function == harmonic_function)
            {
                fprintf(pathToFile, "%.5g \n", sin(position));
            }
            else
            {
                fprintf(pathToFile, "\n");
            }
        }

        double finalStep;
        double nextStep = step;

        if (position + nextStep > rightEndPoint) //Account for interval overshoot
        {
            nextStep = rightEndPoint - position;
        }
        do
        {
            rk_step12(function, position, currentFunctionValue, nextStep, functionStep, functionError);
            error = gsl_blas_dnrm2(functionError);
            normFunctionValue = gsl_blas_dnrm2(functionStep);
            tolerance = (normFunctionValue * relativeAccuracy + absoluteAccuracy) *
                        sqrt(nextStep / (rightEndPoint - leftEndpoint));
            finalStep = nextStep;
            nextStep *= pow(tolerance / error, 0.25) * 0.95; //Gradually make steps smaller
        }
        while (error > tolerance);

        gsl_vector_memcpy(currentFunctionValue, functionStep);
        position += finalStep;

    }

    gsl_vector_memcpy(functionValueRight, functionStep); //Function value at right end point

    //Free memory
    gsl_vector_free(functionStep);
    gsl_vector_free(functionError);
    gsl_vector_free(currentFunctionValue);

}