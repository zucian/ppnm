#include <math.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include "integrationMethods.h"

void print_test_result(char *string, double integralValue, double exactValue, double absoluteAccuracy,
                       double relativeAccuracy, double integrationError, int numberOfCalls)
{
    printf("%s: %g\n", string, integralValue);
    printf("Error goal: %.25g \n", absoluteAccuracy + fabs(exactValue) * relativeAccuracy);
    printf("Actual error: %.25g \n", fabs(integralValue - exactValue));
    printf("Error estimate: %.25g \n", integrationError);
    printf("Function calls: %i \n", numberOfCalls);
}

int main(int argc, char *argv[])
{
    //PART A
    printf("A: Test recursive adaptive integrator: \n\n");

    double leftEndpoint = 0;
    double rightEndpoint = 1;
    double absoluteAccuracy = 1e-3;
    double relativeAccuracy = 1e-3;
    int numberOfCalls = 0;
    double integrationError = 0;

    //First function
    printf("Testing integrator on ∫_0^1 dx √(x) = 2/3 \n\n");
    double exactValue = 2.0 / 3.0;
    double first_test_function(double x)
    {
        numberOfCalls++;
        return sqrt(x);
    }
    double integralValue = integrate(first_test_function, leftEndpoint, rightEndpoint, absoluteAccuracy,
                                     relativeAccuracy, &integrationError);
    print_test_result("Result of numerical integration", integralValue, exactValue, absoluteAccuracy, relativeAccuracy,
                      integrationError, numberOfCalls);


    //Second function
    printf("\nTesting integrator on ∫_0^1 dx 4√(1-x²) = π \n\n");
    exactValue = M_PI;
    numberOfCalls = 0;
    integrationError = 0;
    double second_test_function(double x)
    {
        numberOfCalls++;
        return sqrt(1 - x * x) * 4;
    }
    integralValue = integrate(second_test_function, leftEndpoint, rightEndpoint, absoluteAccuracy, relativeAccuracy,
                              &integrationError);
    print_test_result("Result of numerical integration", integralValue, exactValue, absoluteAccuracy, relativeAccuracy,
                      integrationError, numberOfCalls);
    printf("\n");


    //PART B
    printf("B: Test open quadrature with Clenshaw-Curtis variable transformation\n\n");

    //First function
    printf("Testing Clenshaw-Curtis on ∫_0^1 dx 1/√(x) = 2 \n\n");
    exactValue = 2.0;
    numberOfCalls = 0;
    integrationError = 0;
    double third_test_function(double x)
    {
        numberOfCalls++;
        return 1 / sqrt(x);
    }
    integralValue = open_quad(third_test_function, leftEndpoint, rightEndpoint, absoluteAccuracy,
                              relativeAccuracy, &integrationError);
    print_test_result("Result of numerical integration with Clenshaw-Curtis", integralValue, exactValue,
                      absoluteAccuracy, relativeAccuracy, integrationError, numberOfCalls);
    printf("\n");

    numberOfCalls = 0;
    integrationError = 0;
    integralValue = adapt(third_test_function, leftEndpoint, rightEndpoint, absoluteAccuracy, relativeAccuracy,
                          &integrationError);
    print_test_result("Result of numerical integration without Clenshaw-Curtis", integralValue, exactValue,
                      absoluteAccuracy, relativeAccuracy, integrationError, numberOfCalls);
    printf("\n");

    //Second function
    printf("Testing Clenshaw-Curtis on ∫_0^1 dx ln(x)/√(x) = -4 \n\n");
    exactValue = -4.0;
    numberOfCalls = 0;
    integrationError = 0;
    double fourth_test_function(double x)
    {
        numberOfCalls++;
        return log(x) / sqrt(x);
    }
    integralValue = open_quad(fourth_test_function, leftEndpoint, rightEndpoint, absoluteAccuracy,
                              relativeAccuracy, &integrationError);
    print_test_result("Result of numerical integration with Clenshaw-Curtis", integralValue, exactValue,
                      absoluteAccuracy, relativeAccuracy, integrationError, numberOfCalls);
    printf("\n");

    numberOfCalls = 0;
    integrationError = 0;
    integralValue = adapt(fourth_test_function, leftEndpoint, rightEndpoint, absoluteAccuracy, relativeAccuracy,
                          &integrationError);
    print_test_result("Result of numerical integration without Clenshaw-Curtis", integralValue, exactValue,
                      absoluteAccuracy, relativeAccuracy, integrationError, numberOfCalls);

    printf("\n");

    //Third function
    printf("Testing Clenshaw-Curtis on ∫_0^1 dx 4√(1-x²) = π \n\n");
    exactValue = M_PI;
    numberOfCalls = 0;
    integrationError = 0;
    integralValue = open_quad(second_test_function, leftEndpoint, rightEndpoint, absoluteAccuracy,
                              relativeAccuracy, &integrationError);
    print_test_result("Result of numerical integration with Clenshaw-Curtis", integralValue, exactValue,
                      absoluteAccuracy, relativeAccuracy, integrationError, numberOfCalls);
    printf("\n");

    numberOfCalls = 0;
    integrationError = 0;
    integralValue = adapt(second_test_function, leftEndpoint, rightEndpoint, absoluteAccuracy, relativeAccuracy,
                          &integrationError);
    print_test_result("Result of numerical integration without Clenshaw-Curtis", integralValue, exactValue,
                      absoluteAccuracy, relativeAccuracy, integrationError, numberOfCalls);

    printf("\n");

    double gsl_test_function(double x, void *parameters)
    {
        parameters = NULL;
        return second_test_function(x);
    }

    gsl_function my_gsl_test_function;
    my_gsl_test_function.function = &gsl_test_function;
    my_gsl_test_function.params = NULL;
    size_t limit = 999;
    gsl_integration_cquad_workspace *workspace = gsl_integration_cquad_workspace_alloc(limit);
    double result;
    double absError;
    size_t numOfEvals;

    integrationError = 0;
    gsl_integration_cquad(&my_gsl_test_function, leftEndpoint, rightEndpoint, absoluteAccuracy, relativeAccuracy,
                          workspace,
                          &result, &absError, &numOfEvals);

    print_test_result("Result of numerical integration with GSL Clenshaw-Curtis", result, exactValue, absoluteAccuracy,
                      relativeAccuracy, absError, (int) numOfEvals);
    workspace = NULL;
    gsl_integration_cquad_workspace_free(workspace);

    //double err = 0;

    //PART C
    printf("\nC: Infinite Limit \n\n");

    double fifth_test_function(double x)
    {
        numberOfCalls++;
        return exp(-x * x);
    }

    exactValue = sqrt(M_PI);
    printf("Testing infinite limits on ∫_-inf^inf dx exp(-x²) = √π \n\n");
    numberOfCalls = 0;
    integrationError = 0;
    integralValue = integrate(fifth_test_function, -INFINITY, INFINITY, absoluteAccuracy, relativeAccuracy,
                              &integrationError);
    print_test_result("Result of infinite limit numerical integration", integralValue, exactValue, absoluteAccuracy,
                      relativeAccuracy, integrationError, numberOfCalls);

    printf("\n");

    double gsl_test_function_2(double x, void *params)
    {
        params = NULL;
        return fifth_test_function(x);
    }

    gsl_function my_gsl_test_function_2;
    my_gsl_test_function_2.function = &gsl_test_function_2;
    my_gsl_test_function_2.params = NULL;
    gsl_integration_workspace *workspace2 = gsl_integration_workspace_alloc(limit);
    result = 0;
    absError = 0;
    numOfEvals = 0;
    gsl_integration_qagi(&my_gsl_test_function_2, absoluteAccuracy, relativeAccuracy, limit, workspace2, &result,
                         &absError);
    print_test_result("Result of infinite limit GSL integration", result, exactValue, absoluteAccuracy,
                      relativeAccuracy, absError, (int) numOfEvals);
    printf("QAGI-method has no option to return number of evaluations, hence function calls = 0 \n\n");

    double sixth_test_function(double x)
    {
        numberOfCalls++;
        return 1/(x * x + 1);
    }

    exactValue = M_PI / 2;
    printf("Testing infinite limits on ∫_0^inf dx 1/(1+x²) = π/2 \n\n");
    numberOfCalls = 0;
    integrationError = 0;
    integralValue = integrate(sixth_test_function, 0, INFINITY, absoluteAccuracy, relativeAccuracy,
                              &integrationError);
    print_test_result("Result of infinite limit numerical integration", integralValue, exactValue, absoluteAccuracy,
                      relativeAccuracy, integrationError, numberOfCalls);

    printf("\n");

    double gsl_test_function_3(double x, void *params)
    {
        params = NULL;
        return sixth_test_function(x);
    }

    gsl_function my_gsl_test_function_3;
    my_gsl_test_function_3.function = &gsl_test_function_3;
    my_gsl_test_function_3.params = NULL;
    gsl_integration_workspace *workspace3 = gsl_integration_workspace_alloc(limit);
    result = 0;
    absError = 0;
    numOfEvals = 0;
    gsl_integration_qagiu(&my_gsl_test_function_3, 0, absoluteAccuracy, relativeAccuracy, limit, workspace3, &result,
                          &absError);
    print_test_result("Result of infinite limit GSL integration", result, exactValue, absoluteAccuracy,
                      relativeAccuracy, absError, (int) numOfEvals);
    printf("QAGI-method has no option to return number of evaluations, hence function calls = 0 \n\n");

    return 0;
}
