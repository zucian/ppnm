#ifndef HAVE_MONTECARLO_H
#define HAVE_MONTECARLO_H

void random_numbers(int dimension, const double *lowerBound, const double *upperBound, double *numbers);

double van_der_corput_sequence(int i, int base);

double halton_sequence_first(int i, int dimension, double *vector);

double halton_sequence_second(int i, int dimension, double *vector);

void random_number_halton_corput_first(int i, int dimension, const double *lowerBound, const double *upperBound,
                                       double *vector);

void random_number_halton_corput_second(int i, int dimension, const double *lowerBound, const double *upperBound,
                                        double *vector);

void
plain_monte_carlo(int dimension, double *lowerBound, double *upperBound, double function(double *numbers),
                  int numberOfPoints,
                  double *result, double *error);

void plain_monte_carlo_quasi(int dimension, double *lowerBound, double *upperBound, double function(double *numbers),
                             int numberOfPoints, double *result, double *error);

double plain_monte_carlo_stratified_sampling(int dimension, double function(double *vector), double *lowerBound,
                                             double *upperBound, double absoluteAccuracy, double relativeAccuracy,
                                             int numberOfRecalls, double meanRecalls);

#endif //HAVE_MONTECARLO_H
