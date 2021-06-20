#ifndef HAVE_INTEGRATIONMETHODS_H
#define HAVE_INTEGRATIONMETHODS_H

double adapt_24(double function(double), double leftEndpoint, double rightEndpoint, double absoluteAccuracy,
                double relativeAccuracy, double secondFunctionValue, double thirdFunctionValue, int numberOfRecursions,
                double *IntegrationError);

double adapt(double function(double), double leftEndpoint, double rightEndpoint, double absoluteAccuracy,
             double relativeAccuracy, double *integrationError);

double open_quad(double function(double), double leftEndpoint, double rightEndpoint, double absoluteAccuracy,
                 double relativeAccuracy, double *integrationError);

double integrate(double function(double), double leftEndpoint, double rightEndpoint, double absoluteAccuracy,
                 double relativeAccuracy, double *integrationError);

#endif //HAVE_INTEGRATIONMETHODS_H
