#ifndef HAVE_MINIMIATION_H
#define HAVE_MINIMIATION_H

#include <gsl/gsl_vector.h>

void numerical_gradient(double function(gsl_vector *), gsl_vector *minimum, gsl_vector *gradient);

void quasi_newton_method(double function(gsl_vector *), gsl_vector *minimum, double tolerance);

void simplex_reflection(int dimension, const double *highestPoint, const double *centroidPoint, double *reflectPoint);

void simplex_expansion(int dimension, const double *highestPoint, const double *centroidPoint, double *expandedPoint);

void
simplex_contraction(int dimension, const double *highestPoint, const double *centroidPoint, double *contractedPoint);

void simplex_reduction(int dimension, double **simplex, int lowPointID);

double simplex_distance(int dimension, double *firstPoint, double *secondPoint);

double simplex_size(int dimension, double **simplex);

void simplex_update(int dimension, double **simplex, const double *functionValues, int *highPointID, int lowPointID,
                    double *centroidPoint);

void
simplex_initiate(int dimension, double function(double *), double **simplex, double *functionValues, int *highPointID,
                 int *lowPointID,
                 double *centroidPoint);

int downhill_simplex(int dimension, double function(double *), double **simplex, double simplexSizeGoal);

#endif //HAVE_MINIMIATION_H
