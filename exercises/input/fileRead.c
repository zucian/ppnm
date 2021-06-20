#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "No arguments passed \n");
    }
    else
    {
        int input = 1;
        int X;
        FILE *inputFileStream = fopen(argv[1], "r");
        FILE *outputFileStream = fopen(argv[2], "w");

        while (input != EOF)
        {
            if (input == 1)
            {
                input = fscanf(inputFileStream, "%d", &X);
                fprintf(outputFileStream, "%d %g %g \n", X, sin(X), cos(X));
            }
            else
            {
                fprintf(stderr, "Could not read input \n");
                exit(-1);
            }
        }
        fprintf(stderr, "Program Successful \n");
    }
    return 0;
}