#ifndef HAVE_DIFFERENTIALEQUATIONS_H
#define HAVE_DIFFERENTIALEQUATIONS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void harmonic_function(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative);

void SIR_model(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative);

void SIR_model_new_contact_time(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative);

void three_body_problem(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative);

#endif //HAVE_DIFFERENTIALEQUATIONS_H
