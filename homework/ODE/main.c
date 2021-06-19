#include <math.h>
#include <gsl/gsl_matrix.h>
#include "rungeKutta.h"
#include "differentialEquations.h"

int main(int argc, char *argv[])
{
    //PART A

    //Solving u'' = -u
    int harmonicDimension = 2;
    gsl_vector *harmonicFunctionValueLeft = gsl_vector_alloc(harmonicDimension);
    gsl_vector *harmonicFunctionValueRight = gsl_vector_alloc(harmonicDimension);
    gsl_vector_set(harmonicFunctionValueLeft, 1, 1); //Initial value

    //Settings for ODE solver
    double leftEndpoint = 0.0;
    double rightEndPoint = 3 * M_PI;
    double absoluteAccuracy = 1e-3;
    double relativeAccuracy = 1e-3;
    double step = (rightEndPoint - leftEndpoint) / 10;

    FILE *harmonicOutput = fopen(argv[1], "w");
    rk_driver(&harmonic_function, leftEndpoint, harmonicFunctionValueLeft, rightEndPoint, harmonicFunctionValueRight,
              step, absoluteAccuracy, relativeAccuracy, harmonicOutput);
    fclose(harmonicOutput);

    //Solving SIR-model
    int SIRDimension = 3;
    gsl_vector *SIRFunctionValueLeft = gsl_vector_alloc(SIRDimension);
    gsl_vector *SIRFunctionValueRight = gsl_vector_alloc(SIRDimension);

    leftEndpoint = 0;
    rightEndPoint = 100;

    //COVID-19 data in Denmark 12/4/2021
    double population = 5806000;
    double totalInfected = 237792;
    double recovered = 226630;
    double currentlyInfected = totalInfected - recovered;
    double vaccinated = 445566;
    double dead = 2441;
    double removed = dead + recovered + vaccinated;

    gsl_vector_set(SIRFunctionValueLeft, 0, population - currentlyInfected - removed);
    gsl_vector_set(SIRFunctionValueLeft, 1, currentlyInfected);
    gsl_vector_set(SIRFunctionValueLeft, 2, removed);

    FILE *SIROutput = fopen(argv[2], "w");
    rk_driver(&SIR_model, leftEndpoint, SIRFunctionValueLeft, rightEndPoint, SIRFunctionValueRight, step,
              absoluteAccuracy, relativeAccuracy, SIROutput);
    fclose(SIROutput);

    //SIR again, this time with doubled contact time
    gsl_vector *SIRFunctionValueLeft2 = gsl_vector_alloc(SIRDimension);
    gsl_vector *SIRFunctionValueRight2 = gsl_vector_alloc(SIRDimension);
    gsl_vector_set(SIRFunctionValueLeft2, 0, population - currentlyInfected - removed);
    gsl_vector_set(SIRFunctionValueLeft2, 1, currentlyInfected);
    gsl_vector_set(SIRFunctionValueLeft2, 2, removed);
    FILE *SIROutput2 = fopen(argv[3], "w");
    rk_driver(&SIR_model_new_contact_time, leftEndpoint, SIRFunctionValueLeft2, rightEndPoint, SIRFunctionValueRight2, step,
              absoluteAccuracy, relativeAccuracy, SIROutput2);
    fclose(SIROutput2);

    //PART C

    //Initial values from article https://arxiv.org/abs/math/0011268
    leftEndpoint = 0.0;
    rightEndPoint = 6.32591398;
    int threeBodyDimension = 12;

    gsl_vector *threeBodyFunctionValueLeft = gsl_vector_alloc(threeBodyDimension);
    gsl_vector *threeBodyFunctionValueRight = gsl_vector_alloc(threeBodyDimension);

    double initialX1 = 0.97000436;
    double initialY1 = -0.24308753;
    double initialX2 = -0.97000436;
    double initialY2 = 0.24308753;
    double initialX3 = 0;
    double initialY3 = 0;
    double initialVelocityX3 = -0.93240737;
    double initialVelocityY3 = -0.86473146;
    double initialVelocityX1 = -initialVelocityX3 / 2;
    double initialVelocityY1 = -initialVelocityY3 / 2;
    double initialVelocityX2 = -initialVelocityX3 / 2;
    double initialVelocityY2 = -initialVelocityY3 / 2;

    gsl_vector_set(threeBodyFunctionValueLeft, 0, initialX1);
    gsl_vector_set(threeBodyFunctionValueLeft, 1, initialY1);
    gsl_vector_set(threeBodyFunctionValueLeft, 2, initialX2);
    gsl_vector_set(threeBodyFunctionValueLeft, 3, initialY2);
    gsl_vector_set(threeBodyFunctionValueLeft, 4, initialX3);
    gsl_vector_set(threeBodyFunctionValueLeft, 5, initialY3);
    gsl_vector_set(threeBodyFunctionValueLeft, 6, initialVelocityX1);
    gsl_vector_set(threeBodyFunctionValueLeft, 7, initialVelocityY1);
    gsl_vector_set(threeBodyFunctionValueLeft, 8, initialVelocityX2);
    gsl_vector_set(threeBodyFunctionValueLeft, 9, initialVelocityY2);
    gsl_vector_set(threeBodyFunctionValueLeft, 10, initialVelocityX3);
    gsl_vector_set(threeBodyFunctionValueLeft, 11, initialVelocityY3);

    FILE *threeBodyOutput = fopen(argv[4], "w");
    rk_driver(&three_body_problem, leftEndpoint, threeBodyFunctionValueLeft, rightEndPoint, threeBodyFunctionValueRight,
              step, absoluteAccuracy, relativeAccuracy, threeBodyOutput);
    fclose(threeBodyOutput);

    //Free memory
    gsl_vector_free(harmonicFunctionValueLeft);
    gsl_vector_free(harmonicFunctionValueRight);
    gsl_vector_free(SIRFunctionValueLeft);
    gsl_vector_free(SIRFunctionValueRight);
    gsl_vector_free(SIRFunctionValueLeft2);
    gsl_vector_free(SIRFunctionValueRight2);
    gsl_vector_free(threeBodyFunctionValueLeft);
    gsl_vector_free(threeBodyFunctionValueRight);

    return 0;
}