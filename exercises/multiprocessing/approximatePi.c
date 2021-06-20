#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "randomNumber.h"
#include "approximatePi.h"


#define NUMBER_OF_THREADS 2

typedef struct placePointsStruct
{
    unsigned int seed;
    double numberOfPoints;
    int *pointsInCircle;
} placePointsStruct;

void *place_points(void *placePointsStructInput)
{
    placePointsStruct *arguments = (placePointsStruct *) placePointsStructInput;
    unsigned int seed = arguments->seed;
    double numberOfPoints = arguments->numberOfPoints;
    int *pointsInCircle = arguments->pointsInCircle;
    double unitCircleRadius = 1.0;
    double x = 0;
    double y = 0;

    for (int i = 0; i < numberOfPoints; i++)
    {
        //Generate random points
        x = random_number(&seed);
        y = random_number(&seed);

        if (sqrt(pow(x, 2) + pow(y, 2)) <= unitCircleRadius)
        {
            (*pointsInCircle)++;
        }
    }
    pthread_exit((void *) placePointsStructInput);
    return NULL;
}

void approximate_pi_multi(double numberOfPoints)
{
    pthread_t threadArray[NUMBER_OF_THREADS];
    pthread_attr_t attributes;
    placePointsStruct *placePointsStructArray = malloc(NUMBER_OF_THREADS * sizeof(placePointsStruct));
    int error;
    void *status;

    //Initialize threads
    pthread_attr_init(&attributes);
    pthread_attr_setdetachstate(&attributes, PTHREAD_CREATE_JOINABLE);

    unsigned int masterSeed = time(NULL);
    int *pointsInCircleArray = malloc(NUMBER_OF_THREADS * sizeof(int));

    //Create threads that run in parallel in loop
    for (int i = 0; i < NUMBER_OF_THREADS; i++)
    {
        unsigned int threadSeed = masterSeed + i; //Set unique seed for each thread
        pointsInCircleArray[i] = 0;

        placePointsStructArray[i].seed = threadSeed;
        placePointsStructArray[i].numberOfPoints = numberOfPoints;
        placePointsStructArray[i].pointsInCircle = &pointsInCircleArray[i];

        error = pthread_create(&threadArray[i], &attributes, place_points, (void *) &(placePointsStructArray[i]));
        if (error)
        {
            printf("pthread_join() returned failure");
            exit(-1);
        }

    }
    pthread_attr_destroy(&attributes);

    for (int i = 0; i < NUMBER_OF_THREADS; i++)
    {
        error = pthread_join(threadArray[i], &status);
        if (error)
        {
            printf("pthread_join() returned failure");
            exit(-1);
        }
    }

    double pointsInCircle = 0;
    for (int i = 0; i < NUMBER_OF_THREADS; i++)
    {
        pointsInCircle += pointsInCircleArray[i];
    }
    pointsInCircle /= NUMBER_OF_THREADS;

    double approximatePi = 4.0*pointsInCircle/((double)numberOfPoints);

    printf("Approximation of pi using pthreads = %g \n",approximatePi);
    printf("Deviation from actual pi: %g percent \n", fabs(1.0-approximatePi/M_PI)*100);

    free(pointsInCircleArray);
}