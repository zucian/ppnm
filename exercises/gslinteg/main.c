#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double Afunction(double x, void* params){
    double value = log(x)/(sqrt(x));
    return value;
}

double aIntegral(){
    gsl_function F;
    F.function = &Afunction;
    int limit = 999;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
    double a=0, b=1, acc=1e-6, eps=1e-6, result, error;
    gsl_integration_qags(&F, a, b, acc, eps, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}

int main(){
    FILE* outA = fopen("out.valueExerciseA.txt","w");
    double exerciseA = aIntegral();
    fprintf(outA, "int(log(x)/sqrt(x), x=0..1) = %g\n", exerciseA);
    printf("Result can be seen in out.valueExerciseA.txt file in directory\n");
    return 0;
}