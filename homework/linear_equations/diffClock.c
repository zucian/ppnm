//
// Implemented by Marc B. Sørensen on 2/13/21.
//

#include <time.h>
#include "diffClock.h"
/*
  USE LIKE THIS:
  // -----------------------------------------------------------------------------------------------------------

  clock_t begin = clock(); // We define variables to hold the times used for timing the computations.
  clock_t end   = clock(); // These are defined in <time.h> and used in diffClock(). It is my own implementation.

  begin = clock(); // Begin timing

  // time something....

  end = clock(); // end timing
  printf("Done! Elapsed time = %g ms \n", (double)(diffClock(end, begin)) );


*/

double diffClock(clock_t startTime, clock_t endTime){
  /*
  Function diffClock, takes two arguments, the time stamp
  of the start and end of a time interval of interest.

  ¤ startTime : Timestamp of beginning of time interval
  ¤ endTime   : Timestamp of ending of time interval

  */

  double diffTicks  =  startTime - endTime;
  double diffms     =  (diffTicks * 10) / CLOCKS_PER_SEC;

  return diffms;
}
