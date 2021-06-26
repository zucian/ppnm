#include <stdio.h>
#include <stdlib.h>

#include "utilities.h"

int binarySearch( int numOfPts, double* pts, double evalPt ){
    //  Define indices of interval bounds
    int leftId   =  0 ;
    int rightId  =  numOfPts - 1 ;

    while ( rightId - leftId > 1 ) {        // While there are two intervals to choose from
        int middleId  =  (leftId + rightId) / 2;

        if ( evalPt > pts[middleId] ) {       // If the point is in the right interval
            leftId   =  middleId;
        }
        else {                                // Else it is in the left interval
            rightId  =  middleId;
        }
    }
    int whichInterval = leftId;

    return whichInterval;
}


void inputToArray(double* XData, double* YData, char* inputFilename ){
    int numOfReturnVals = 2;  // The return value of fscanf() is the number of fetched variables

    int    input = numOfReturnVals; 			// Input will receive the return value from fscanf
    double argX;						              // Variable to hold the X data from stdin
    double argY;						              // Variable to hold the X data from stdin

    FILE* myInputFileStream  = fopen(inputFilename, "r");

    int id = 0;
    while( input != EOF) { // While we are not at the end of the input stream
        if ( input == numOfReturnVals){		 // If input has returned success
            input = fscanf(myInputFileStream, "%lg\t%lg", &argX, &argY); // Fetch data from stream and place it at the address of argX

            XData[id] = argX;
            YData[id] = argY;

            id++;
        }
        else { // If we are not successfull, say if there was an error with fscanf()
            fprintf(stderr, "Failed to read input.\n"); // Print this to stderr so we know
            exit(-1); 																	// Terminate program
        }
    }

    fclose(myInputFileStream);
}

void inputToArray3D(double* XData, double* YData, double* ZData, char* inputFilename ){
    int numOfReturnVals = 2;  // The return value of fscanf() is the number of fetched variables

    int    input = numOfReturnVals; 			// Input will receive the return value from fscanf
    double argX;						              // Variable to hold the X data from stdin
    double argY;						              // Variable to hold the X data from stdin
    double argZ;

    FILE* myInputFileStream  = fopen(inputFilename, "r");

    int id = 0;
    while( input != EOF) { // While we are not at the end of the input stream
        if ( input == numOfReturnVals){		 // If input has returned success
            input = fscanf(myInputFileStream, "%lg\t%lg\t%lg", &argX, &argY,&argZ); // Fetch data from stream and place it at the address of argX

            XData[id] = argX;
            YData[id] = argY;
            ZData[id] = argZ;

            id++;
        }
        else { // If we are not successfull, say if there was an error with fscanf()
            fprintf(stderr, "Failed to read input.\n"); // Print this to stderr so we know
            exit(-1); 																	// Terminate program
        }
    }

    fclose(myInputFileStream);
}