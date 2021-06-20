#include <stdlib.h>
#include "randomNumber.h"

double random_number(unsigned int *seed)
{
    double maxRandom = (double) RAND_MAX;
    double randomNumber = (double) rand_r(seed);

    return 2*randomNumber/maxRandom-1; //rescale between -1 and 1
}