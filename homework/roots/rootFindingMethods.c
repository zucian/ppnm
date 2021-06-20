#include <float.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "rootFindingMethods.h"
#include "GS_utilities.h"

void
newton_method(void function(gsl_vector *point, gsl_vector *functionValues), gsl_vector *startingPoint, double tolerance)
{
    int count = 0;
    int dimension = startingPoint->size;
    double stepSize = sqrt(DBL_EPSILON);

    //Allocate memory
    gsl_matrix *jacobianMatrix = gsl_matrix_alloc(dimension, dimension);
    gsl_matrix *triangleMatrix = gsl_matrix_alloc(dimension, dimension);
    gsl_vector *functionValue = gsl_vector_alloc(dimension);
    gsl_vector *functionValueTemporary = gsl_vector_alloc(dimension);
    gsl_vector *solution = gsl_vector_alloc(dimension);
    gsl_vector *solutionScaled = gsl_vector_alloc(dimension);
    gsl_vector *nextFunctionValue = gsl_vector_alloc(dimension);
    gsl_vector *nextPoint = gsl_vector_alloc(dimension);

    function(startingPoint, functionValue);
    while (gsl_blas_dnrm2(functionValue) > tolerance)
    {
        count++;
        assert(count < 1e5);
        for (int i = 0; i < dimension; i++)
        {
            gsl_vector_set(startingPoint, i, gsl_vector_get(startingPoint, i) + stepSize);
            function(startingPoint, functionValueTemporary);
            for (int j = 0; j < dimension; j++)
            {
                double functionValueTemporaryI = gsl_vector_get(functionValueTemporary, j);
                double functionValueI = gsl_vector_get(functionValue, j);
                double functionValueDifference = functionValueTemporaryI - functionValueI;
                double matrixElement = functionValueDifference / stepSize;
                gsl_matrix_set(jacobianMatrix, j, i, matrixElement);
            }
            gsl_vector_set(startingPoint, i, gsl_vector_get(startingPoint, i) - stepSize);
        }
        gsl_vector_scale(functionValue, -1.0);

        GS_decomp(jacobianMatrix, triangleMatrix);
        GS_solve(jacobianMatrix, triangleMatrix, functionValue, solution);

        gsl_vector_scale(functionValue, -1.0);

        double scale = 2;
        while ((gsl_blas_dnrm2(nextFunctionValue) >= (1 - scale / 2) * gsl_blas_dnrm2(functionValue)) && scale >= 0.02)
        {
            scale /= 2;

            gsl_vector_memcpy(solutionScaled, solution);
            gsl_vector_scale(solutionScaled, scale);
            gsl_vector_memcpy(nextPoint, solutionScaled);
            gsl_vector_add(nextPoint, startingPoint);
            function(nextPoint, nextFunctionValue);
        }

        gsl_vector_memcpy(startingPoint,nextPoint);
        gsl_vector_memcpy(functionValue, nextFunctionValue);

        if(gsl_blas_dnrm2(solution) <stepSize){
            break;
        }
    }

    //Free memory
    gsl_matrix_free(jacobianMatrix);
    gsl_matrix_free(triangleMatrix);
    gsl_vector_free(functionValue);
    gsl_vector_free(functionValueTemporary);
    gsl_vector_free(solution);
    gsl_vector_free(solutionScaled);
    gsl_vector_free(nextFunctionValue);
    gsl_vector_free(nextPoint);
}