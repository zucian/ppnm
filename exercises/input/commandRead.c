#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    if (argc < 2) //Check for passed arguments
    {
        fprintf(stderr, "No arguments passed \n");
    }
    else
    {
        for (int i = 1; i < argc; i++)
        {
            double X = atof(argv[i]);
            fprintf(stdout, "%g %g %g \n", X, sin(X), cos(X));
        }
        fprintf(stderr, "Program Successful \n");
    }
    return 0;
}