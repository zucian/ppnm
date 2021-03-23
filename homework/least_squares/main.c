#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

typedef double (*f)(double x);

double f_0(double x){
    return 1.;
}

double f_1(double x){
    return x;
}

void least_sqr_fit(gsl_matrix* dataMatrix, f functions[],int nFunctions, gsl_vector* c){
    gsl_matrix* A = gsl_matrix_alloc(dataMatrix->size1,nFunctions);
    gsl_matrix* R = gsl_matrix_alloc(dataMatrix->size1,nFunctions);
    gsl_vector* y = gsl_vector_alloc(dataMatrix->size1);
    gsl_matrix_get_col(y, dataMatrix, 1);

    for(int i=0; i<(A->size1); i++){
        double x_i = gsl_matrix_get(data_matrix,i,0); //Gets the i'th x-value and stores it in x_i
        for (int j=0; j<(A->size2); j++){
            double A_ij = functions[j](x_i); //Computes f_j(x_i) and stores it in A_ij
            gsl_matrix_set(A,i,j,A_ij); //Set (i,j) in A equal to A_ij
        }
    }
    GS_decomp(A,R); //Gram-schmidt storing orthogonal Q in A and R is upper triangular
    GS_solve(A,R,y,c) //Solves R*c = Q^t*b and stores in c

    gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_vector_free(y);
}

int main(){
    int nData = 9;
    int nFunctions = 2;
    double x[nData] = {1., 2., 3., 4., 6., 9., 10., 13., 15.};
    double y[nData] = {117., 100., 88., 72.,  53., 29.5, 25.2, 15.2, 11.1};

    gsl_matrix* dataMatrix = gsl_matrix_alloc(nData,3);
    gsl_vector* c = gsl_vector_alloc(nData);
    for(int i=0; i<nDatapoints; i++){
        gsl_matrix_set(dataMatrix,i,0,x[i]); //Set first column in dataMatrix equal to x-values
        gsl_matrix_set(dataMatrix,i,1,log(y[i])); //Set second column in dataMatrix equal to log of y-values
        gsl_matrix_set(dataMatrix,i,2,(y[i]/20)/y[i]); //Set third column in dataMatrix equal to uncertainty in y
    }

    functions[2] = {&f0, &f1};

    least_sqr_fit(dataMatrix,functions,nFunctions,c);

    return 0;
}