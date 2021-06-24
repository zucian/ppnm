#include <gsl/gsl_integration.h>
#include <math.h>

#include "integrateFunction.h"

double integrate_function( double lowerLimit, double upperLimit, const gsl_function* myGSLFUNC, double tolAbs, double tolRel, size_t limit ){
    /* Function to integrate the GSL_FUNCTION on (0, 1)
     *
     *  ¤ double lowerLimit             : lower integration limit.
     *  ¤ double upperLimit             : upper integration limit.
     *  ¤ const gsl_function* myGSLFUNC : A gsl_function* struct
     *  ¤ double tolAbs                 : Tolerance for absolute error estimate
     *  ¤ double tolRel                 : Tolerance for absolute error estimate
     *  ¤ size_t limit                  : The maximum number of iterations
     *  ¤
     */

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc( limit );
    double result;
    double absError;

    gsl_integration_qags(myGSLFUNC, lowerLimit, upperLimit, tolAbs, tolRel, limit, workspace, &result, &absError);

    gsl_integration_workspace_free(workspace);
    workspace = NULL;

    return result;
}