#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv){
    if(argc<2) fprintf(stderr, "No arguments were passed \n");
    else{
        for (int i=1; i<argc; i++){
            double x = atof(argv[i]); // convert string to floating point
            fprintf(stdout,"%g %g %g\n",x,sin(x),cos(x)); // print to stdout
        }
    }
    return 0;
}