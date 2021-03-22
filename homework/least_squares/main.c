#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

void least_sqr_fit(gsl_matrix* data_matrix, f functions[], gsl_vector* c){
    gsl_matrix* A = gsl_matrix_alloc(data_matrix->size1,function_output);
    gsl_vector* y = gsl_vector_alloc(data_matrix->size1);
    gsl_matrix_get_col(y, data_matrix, 1);

    for(int i=0; i<(A->size1); i++){
        double x_i = gsl_matrix_get(data_matrix,i,0); //Gets the i'th x-value and stores it in x_i
        for (int j=0; j<(A->size2); j++){
            double A_ij = functions[j](x_i); //Computes f_j(x_i) and stores it in A_ij
            gsl_matrix_set(A,i,j,A_ij);
        }
    }
    GS_decomp(A,R); //Gram-schmidt storing orthogonal Q in A and R is upper triangular
    GS_solve(A,R,y,c) //Solves R*c = Q^t*b and stores in c

    gsl_matrix_free(A);
    gsl_matrix_free(R);
}

int main(){
    double x[] = {1., 2., 3., 4., 6., 9., 10., 13., 15.};
    double y[] = {117., 100., 88., 72.,  53., 29.5, 25.2, 15.2, 11.1};

    return 0;
}