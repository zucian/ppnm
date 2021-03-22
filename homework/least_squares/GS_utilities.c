#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

void matrix_print(gsl_matrix* A,FILE* file){
    for(int i=0; i<(A->size1);i++){
        for(int j=0; j<(A->size2); j++){
            double Aij = gsl_matrix_get(A,i,j);
            fprintf(file, "%0.3g ",Aij);
        }
        fprintf(file,"\n");
    }
}

void GS_decomp(gsl_matrix* A, gsl_matrix* R){

    //Takes a matrix A and performs Gram-Schmidt orthogonalization. On exit A is replaced with orthogonal matrix Q and also returns upper triangular matrix R

    //size1 refers to amount of rows. size2 refers to columns.
    gsl_vector* a_i = gsl_vector_alloc((A->size1)); //Allocate memory for a_i and q_i to be used later.
    gsl_vector* q_i = gsl_vector_alloc((A->size1)); //Size should be equal to amount of rows in A

    for (int i=0; i<(A->size2); i++){ //Loop runs over every column of A, hence size2
        gsl_matrix_get_col(a_i, A, i); //Copies the i'th column of A into a_i
        gsl_vector_memcpy(q_i,a_i); //Copies a_i into q_i

        double a_i_norm = gsl_blas_dnrm2(a_i); //Computes the norm of a_i
        gsl_matrix_set(R,i,i, a_i_norm); //Set R_ii to norm of a_i
        gsl_vector_scale(q_i,1/a_i_norm); //Scale q_i with 1/||a_i|| so q_i = a_i/||a_i||

        gsl_matrix_set_col(A,i,q_i); //Set i'th column of A equal to q_i

        for (int j=i+1; j<(A->size2); j++){ //Loop runs over columns to the right of i
            gsl_matrix_get_col(a_i,A,j); //Copies the i'th column of A into a_j. a_i and a_j are the same in this case

            double R_ij=0; //Define double for the (i,j) entrance of matrix R
            gsl_blas_ddot(q_i,a_i,&R_ij); //Evaluates is q_i^t*a_i and stores the value in R_ij
            gsl_matrix_set(R,i,j,R_ij); //Set (i,j) entrance of R to R_ij
            gsl_matrix_set(R,j,i,0.); //Set the (j,i) entrance of R to 0. These are the lower triangular elements of R that we want to be zero.

            gsl_blas_daxpy(-R_ij,q_i,a_i); //Returns -R_ij*q_i+a_i and stores it in a_i
            gsl_matrix_set_col(A,j,a_i); //Set j'th column of A equal to a_i
        }
    }
    gsl_vector_free(a_i);
    gsl_vector_free(q_i);
}

void backsub(gsl_matrix* R, gsl_vector* x){

    //Takes a upper triangular matrix R and performs backsubstitution. Returns solution in x.

    for(int i=(x->size)-1; i>=0; i--){ //Loop runs from last column of R backwards
        double R_ii = gsl_matrix_get(R,i,i); //Stores the i'th diagonal entrance of R in R_ii
        double x_i = gsl_vector_get(x,i); //Stores the i'th entrance of x in x_i

        for (int j=i+1; j<(x->size);j++){ //Loop runs over rows above i'th diagonal entrance
            double R_ij = gsl_matrix_get(R,i,j); //Stores (i,j) entrance of R in R_ij
            double x_j = gsl_vector_get(x,j); //Stores the j'th entrance of x in x_j
            x_i -= R_ij*x_j;
        }
        x_i /= R_ii;
        gsl_vector_set(x,i,x_i); //Set the i'th entrance of x equal to x_i
    }
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){

    //Solves the equation  QRx=b by applying Q^T to the vector b and then performing back_substitution on x

    gsl_blas_dgemv(CblasTrans, 1.0, Q,b,0.0,x); //Computes the matrix vector product Q^T*b and saves the result in x
    backsub(R,x); //Compute backsubstitution of R and x and returns answer in x
}