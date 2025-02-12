#include"stdio.h"
#include<gsl/gsl_matrix.h> //wrong way to include

int print_half_00(gsl_matrix *m)
{
    double half = 1.0 / 2; //Cant divide ints
    int status = printf("half m_{00} = %g\n", gsl_matrix_get(m, 0, 0) * half); //use %g for double
    return status;
}

int main(void)
{
    gsl_matrix *m = gsl_matrix_alloc(1, 1); //cant allocate 0,0
    gsl_matrix_set(m, 0, 0, 66);
    printf("half m_{00} (should be 33):\n");
    int status = print_half_00(m);
    if (status < 0) //smaller than not greater than
        printf("status=%d : SOMETHING WENT TERRIBLY WRONG (status<0)\n", status); //%d for int
    else
        printf("status=%d : everything went just fine (status>=0)\n", status);
    gsl_matrix_free(m);
    return 0;
}
