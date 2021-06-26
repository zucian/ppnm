#ifndef HAVE_UTILITIES_H
#define HAVE_UTILITIES_H

#include <stdio.h>
#include <stdlib.h>

int binarySearch( int numOfPts, double* pts, double evalPt );
void inputToArray(double* XData, double* YData, char* inputFilename );
void inputToArray3D(double* XData, double* YData, double* ZData, char* inputFilename );

#endif