#ifndef UNI_APPROXIMATEPI_H
#define UNI_APPROXIMATEPI_H

#include <pthread.h>

void *place_points(void *placePointsStructInput);

void approximate_pi_multi(double numberOfPoints);

#endif //UNI_APPROXIMATEPI_H
