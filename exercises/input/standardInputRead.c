#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    int input = 1;
    int X;

    while (input != EOF)
    {
        if (input == 1)
        {
            input = fscanf(stdin, "%d", &X);
            fprintf(stdout, "%d %g %g \n", X, sin(X), cos(X));
        }
        else
        {
            fprintf(stderr, "Could not read input \n");
            exit(-1);
        }
    }
    fprintf(stderr, "Program successful \n");

    return 0;
}