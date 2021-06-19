#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<float.h>
#include "minimization.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define MY_DOUBLE_EPSILON 2.22045e-10


void numerical_gradient(double function(gsl_vector *), gsl_vector *minimum, gsl_vector *gradient)
{
    long double stepSize = MY_DOUBLE_EPSILON; //Long double?
    double functionValue = function(minimum);
    int dimension = minimum->size;

    for (int i = 0; i < dimension; i++)
    {
        double step;
        double minimumI = gsl_vector_get(minimum, i);

        if (fabs(minimumI) < stepSize)
        {
            step = stepSize;
        }
        else
        {
            step = fabs(minimumI) * stepSize;
        }

        gsl_vector_set(minimum, i, minimumI + step);
        gsl_vector_set(gradient, i, (function(minimum) - functionValue) / step);
        gsl_vector_set(minimum, i, minimumI - step);
    }
}

void quasi_newton_method(double function(gsl_vector *), gsl_vector *minimum, double tolerance)
{
    double stepSize = MY_DOUBLE_EPSILON;
    int dimension = minimum->size;
    int numberOfSteps = 0;
    int numberOfScales = 0;
    int numberOfResets = 0;

    gsl_matrix *inverseHessianMatrix = gsl_matrix_alloc(dimension, dimension);
    gsl_matrix *identityMatrix = gsl_matrix_alloc(dimension, dimension);
    gsl_matrix_set_identity(inverseHessianMatrix);
    gsl_matrix_set_identity(identityMatrix);

    //Allocate memory for needed parts
    gsl_vector *gradientValue = gsl_vector_alloc(dimension);
    gsl_vector *nextGradientValue = gsl_vector_alloc(dimension);
    gsl_vector *newtonStep = gsl_vector_alloc(dimension);
    gsl_vector *nextMinimum = gsl_vector_alloc(dimension);
    gsl_vector *solution = gsl_vector_alloc(dimension);
    gsl_vector *solutionChange = gsl_vector_alloc(dimension);
    gsl_vector *broydenVector = gsl_vector_alloc(dimension);

    numerical_gradient(function, minimum, gradientValue);
    double functionValue = function(minimum);
    double nextFunctionValue;

    while (numberOfSteps < 1e4)
    {
        numberOfSteps++;
        gsl_blas_dgemv(CblasNoTrans, -1, inverseHessianMatrix, gradientValue, 0, newtonStep);
        if (gsl_blas_dnrm2(newtonStep) < stepSize * gsl_blas_dnrm2(minimum))
        {
            fprintf(stderr, "Quasi_newton_method: |dx| < stepSize*|x|\n");
            break;
        }
        if (gsl_blas_dnrm2(gradientValue) < tolerance)
        {
            fprintf(stderr, "Quasi_newton_method: |grad| < accuracy\n");
            break;
        }

        double scale = 1;

        while (1) //while(1) means running until explict break
        {
            gsl_vector_memcpy(nextMinimum, minimum);
            gsl_vector_add(nextMinimum, newtonStep);
            nextFunctionValue = function(nextMinimum);
            double symmetricTransGradient;
            gsl_blas_ddot(newtonStep, gradientValue, &symmetricTransGradient);

            if (nextFunctionValue < functionValue + 0.01 * symmetricTransGradient)
            {
                numberOfScales++;
                break;
            }
            if (scale < stepSize)
            {
                numberOfResets++;
                gsl_matrix_set_identity(inverseHessianMatrix);
                break;
            }
            scale *= 0.5;
            gsl_vector_scale(newtonStep, 0.5);
        }

        numerical_gradient(function, nextMinimum, nextGradientValue);
        gsl_vector_memcpy(solution, nextGradientValue);
        gsl_blas_daxpy(-1, gradientValue, solution);
        gsl_vector_memcpy(solutionChange, newtonStep);
        gsl_blas_dgemv(CblasNoTrans, -1, inverseHessianMatrix, solution, 1, solutionChange);

        gsl_matrix *solutionChangeSolutionChangeTrans = gsl_matrix_calloc(dimension, dimension); //u*u^t
        gsl_blas_dsyr(CblasUpper, 1.0, solutionChange, solutionChangeSolutionChangeTrans);
        double solutionChangeTransSolution; //u^T*y
        gsl_blas_ddot(solutionChange, solution, &solutionChangeTransSolution);
        if (fabs(solutionChangeTransSolution) > 1e-12)
        {
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0 / solutionChangeTransSolution,
                           solutionChangeSolutionChangeTrans, identityMatrix, 1.0, inverseHessianMatrix);
        }

        gsl_vector_memcpy(minimum, nextMinimum);
        gsl_vector_memcpy(gradientValue, nextGradientValue);
        functionValue = nextFunctionValue;
    }

    //Free allocated memory
    gsl_matrix_free(inverseHessianMatrix);
    gsl_matrix_free(identityMatrix);
    gsl_vector_free(gradientValue);
    gsl_vector_free(nextGradientValue);
    gsl_vector_free(newtonStep);
    gsl_vector_free(nextMinimum);
    gsl_vector_free(solution);
    gsl_vector_free(solutionChange);
    gsl_vector_free(broydenVector);

    fprintf(stderr,
            "Quasi_newton_method: \n amount of steps = %i \n amount of scales = %i, \n amount of matrix resets = %i \n  f(x) = %.1e\n\n",
            numberOfSteps, numberOfScales, numberOfResets, functionValue);
}

void simplex_reflection(int dimension, const double *highestPoint, const double *centroidPoint, double *reflectPoint)
{
    for (int i = 0; i < dimension; i++)
    {
        reflectPoint[i] = 2 * centroidPoint[i] - highestPoint[i];
    }
}

void simplex_expansion(int dimension, const double *highestPoint, const double *centroidPoint, double *expandedPoint)
{
    for (int i = 0; i < dimension; i++)
    {
        expandedPoint[i] = 3 * centroidPoint[i] - 2 * highestPoint[i];
    }
}

void
simplex_contraction(int dimension, const double *highestPoint, const double *centroidPoint, double *contractedPoint)
{
    for (int i = 0; i < dimension; i++)
    {
        contractedPoint[i] = 0.5 * centroidPoint[i] + 0.5 * highestPoint[i];
    }
}

void simplex_reduction(int dimension, double **simplex, int lowPointID)
{
    for (int i = 0; i < dimension + 1; i++)
    {
        if (i != lowPointID)
        {
            for (int j = 0; j < dimension; j++)
            {
                simplex[i][j] = 0.5 * (simplex[i][j] + simplex[lowPointID][j]);
            }
        }
    }
}

double simplex_distance(int dimension, double *firstPoint, double *secondPoint)
{
    double euclideanDistance = 0;
    for (int i = 0; i < dimension; i++)
    {
        euclideanDistance += pow(secondPoint[i] - firstPoint[i], 2);
    }
    return sqrt(euclideanDistance);
}

double simplex_size(int dimension, double **simplex)
{
    double temporaryDistance = 0;
    for (int i = 1; i < dimension + 1; i++)
    {
        double distance = simplex_distance(dimension, simplex[0], simplex[i]);
        if (distance < temporaryDistance)
        {
            temporaryDistance = distance;
        }
    }
    return temporaryDistance;
}

void simplex_update(int dimension, double **simplex, const double *functionValues, int *highPointID, int *lowPointID,
                    double *centroidPoint)
{
    *highPointID = 0;
    *lowPointID = 0;

    double lowestPoint = functionValues[0];
    double highestPoint = functionValues[0];

    for (int i = 1; i < dimension + 1; i++)// Find high/low point
    {
        double nextPoint = functionValues[i];
        if (nextPoint > highestPoint)
        {
            highestPoint = nextPoint;
            *highPointID = i;
        }
        if (nextPoint < lowestPoint)
        {
            lowestPoint = nextPoint;
            *lowPointID = i;
        }
    }

    for (int i = 0; i < dimension; i++)//Find center point
    {
        double sum = 0;
        for (int j = 0; j < dimension + 1; j++)
        {
            if (j != *highPointID)
            {
                sum += simplex[j][i];
            }
        }
        centroidPoint[i] = sum / dimension;
    }
}

void
simplex_initiate(int dimension, double function(double *), double **simplex, double *functionValues, int *highPointID,
                 int *lowPointID,
                 double *centroidPoint)
{
    for (int i = 0; i < dimension + 1; i++)
    {
        functionValues[i] = function(simplex[i]);
    }
    simplex_update(dimension, simplex, functionValues, highPointID, lowPointID, centroidPoint);
}

int downhill_simplex(int dimension, double function(double *), double **simplex, double simplexSizeGoal)
{
    double centroid[dimension];
    double functionValues[dimension + 1];
    double reflectedPoint[dimension];
    double expandedPoint[dimension];
    int numberOfSteps = 0;
    int highPointID;
    int lowPointID;

    simplex_initiate(dimension, function, simplex, functionValues, &highPointID, &lowPointID, centroid);

    while (simplex_size(dimension, simplex) > simplexSizeGoal)
    {
        simplex_update(dimension, simplex, functionValues, &highPointID, &lowPointID, centroid);

        simplex_reflection(dimension, simplex[highPointID], centroid, reflectedPoint); //Start with reflection
        double functionValueReflected = function(reflectedPoint);
        if (functionValueReflected < functionValues[lowPointID]) //Good reflection followed by expansion
        {
            simplex_expansion(dimension, simplex[highPointID], centroid, expandedPoint);
            double functionValueExpanded = function(expandedPoint);

            if (functionValueExpanded < functionValueReflected) //Good expansion
            {
                for (int i = 0; i < dimension; i++)
                {
                    simplex[highPointID][i] = expandedPoint[i];
                }
                functionValues[highPointID] = functionValueExpanded;
            }
            else // Bad expansion, accept reflection
            {
                for (int i = 0; i < dimension; i++)
                {
                    simplex[highPointID][i] = reflectedPoint[i];
                }
                functionValues[highPointID] = functionValueReflected;
            }
        }
        else //Bad reflection
        {
            if (functionValueReflected < functionValues[highPointID])
            {
                for (int i = 0; i < dimension; i++)
                {
                    simplex[highPointID][i] = reflectedPoint[i];
                }
                functionValues[highPointID] = functionValueReflected;
            }
            else //Contraction
            {
                simplex_contraction(dimension, simplex[highPointID], centroid, reflectedPoint);
                double functionValueContracted = function(reflectedPoint);

                if (functionValueContracted < functionValues[highPointID]) //Good contraction
                {
                    for (int i = 0; i < dimension; i++)
                    {
                        simplex[highPointID][i] = reflectedPoint[i];
                    }
                    functionValues[highPointID] = functionValueContracted;
                }
                else //Do reduction
                {
                    simplex_reduction(dimension, simplex, lowPointID);
                    simplex_initiate(dimension, function, simplex, functionValues, &highPointID, &lowPointID, centroid);
                }
            }
        }
        numberOfSteps++;
    }
    return numberOfSteps;
}

