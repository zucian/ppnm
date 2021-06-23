#ifndef HAVE_UTILITIES_H
#define HAVE_UTILITIES_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

double random_number(unsigned int *seed);

void vector_print(char *string, gsl_vector *vector);

void print_matrix(int numOfRows, gsl_matrix *matrixToPrint, char *string);

void set_data_symmetric(gsl_matrix *inputMatrix, unsigned int *seed);

#endif
