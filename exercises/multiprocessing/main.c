#include <pthread.h>
#include <stdio.h>
#include "approximatePi.h"

int main()
{
    double numberOfPoints = 1e7;

        approximate_pi_multi(numberOfPoints);

    return 0;
}