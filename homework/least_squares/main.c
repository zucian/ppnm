#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

typedef double (*f)(double x);

double f_0(double x)
{
    return 1.;
}

double f_1(double x)
{
    return x;
}

void matrix_print(gsl_matrix *A, FILE *file)
{
    for (int i = 0; i < (A->size1); i++)
    {
        for (int j = 0; j < (A->size2); j++)
        {
            double Aij = gsl_matrix_get(A, i, j);
            fprintf(file, "%0.3g ", Aij);
        }
        fprintf(file, "\n");
    }
}

void GS_decomp(gsl_matrix *A, gsl_matrix *R)
{

    //Takes a matrix A and performs Gram-Schmidt orthogonalization. On exit A is replaced with orthogonal matrix Q and also returns upper triangular matrix R

    //size1 refers to amount of rows. size2 refers to columns.
    gsl_vector *a_i = gsl_vector_alloc((A->size1)); //Allocate memory for a_i and q_i to be used later.
    gsl_vector *q_i = gsl_vector_alloc((A->size1)); //Size should be equal to amount of rows in A

    for (int i = 0; i < (A->size2); i++)
    { //Loop runs over every column of A, hence size2
        gsl_matrix_get_col(a_i, A, i); //Copies the i'th column of A into a_i
        gsl_vector_memcpy(q_i, a_i); //Copies a_i into q_i

        double a_i_norm = gsl_blas_dnrm2(a_i); //Computes the norm of a_i
        gsl_matrix_set(R, i, i, a_i_norm); //Set R_ii to norm of a_i
        gsl_vector_scale(q_i, 1 / a_i_norm); //Scale q_i with 1/||a_i|| so q_i = a_i/||a_i||

        gsl_matrix_set_col(A, i, q_i); //Set i'th column of A equal to q_i

        for (int j = i + 1; j < (A->size2); j++)
        { //Loop runs over columns to the right of i
            gsl_matrix_get_col(a_i, A, j); //Copies the i'th column of A into a_j. a_i and a_j are the same in this case

            double R_ij = 0; //Define double for the (i,j) entrance of matrix R
            gsl_blas_ddot(q_i, a_i, &R_ij); //Evaluates is q_i^t*a_i and stores the value in R_ij
            gsl_matrix_set(R, i, j, R_ij); //Set (i,j) entrance of R to R_ij
            gsl_matrix_set(R, j, i,
                           0.); //Set the (j,i) entrance of R to 0. These are the lower triangular elements of R that we want to be zero.

            gsl_blas_daxpy(-R_ij, q_i, a_i); //Returns -R_ij*q_i+a_i and stores it in a_i
            gsl_matrix_set_col(A, j, a_i); //Set j'th column of A equal to a_i
        }
    }
    gsl_vector_free(a_i);
    gsl_vector_free(q_i);
}

void backsub(gsl_matrix *R, gsl_vector *x)
{

    //Takes a upper triangular matrix R and performs backsubstitution. Returns solution in x.

    for (int i = (x->size) - 1; i >= 0; i--)
    { //Loop runs from last column of R backwards
        double R_ii = gsl_matrix_get(R, i, i); //Stores the i'th diagonal entrance of R in R_ii
        double x_i = gsl_vector_get(x, i); //Stores the i'th entrance of x in x_i

        for (int j = i + 1; j < (x->size); j++)
        { //Loop runs over rows above i'th diagonal entrance
            double R_ij = gsl_matrix_get(R, i, j); //Stores (i,j) entrance of R in R_ij
            double x_j = gsl_vector_get(x, j); //Stores the j'th entrance of x in x_j
            x_i -= R_ij * x_j;
        }
        x_i /= R_ii;
        gsl_vector_set(x, i, x_i); //Set the i'th entrance of x equal to x_i
    }
}

void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x)
{

    //Solves the equation  QRx=b by applying Q^T to the vector b and then performing back_substitution on x, returning answer in x

    gsl_blas_dgemv(CblasTrans, 1.0, Q, b, 0.0, x); //Computes the matrix vector product Q^T*b and saves the result in x
    backsub(R, x); //Compute backsubstitution of R and x and returns answer in x
}

void GS_inverse(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B)
{
    gsl_vector *unitVector = gsl_vector_alloc(R->size1);
    gsl_matrix *RInverse = gsl_matrix_alloc(R->size1, R->size2);
    for (int i = 0; i < (R->size2); i++)
    { //Loops over every column in R
        gsl_vector_set_basis(unitVector, i);
        backsub(R, unitVector); //Solve R*x = e_i, where e_i is the i'th standard basis vector
        gsl_matrix_set_col(RInverse, i, unitVector); //Turns i'th column in R^-1 into solution to R*x = e_i
    }
    gsl_vector_free(unitVector);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., RInverse, Q, 0., B);
}

void least_sqr_fit(gsl_matrix *dataMatrix, f functions[], int nFunctions, gsl_vector *c, gsl_matrix *covMatrix)
{
    gsl_matrix *A = gsl_matrix_alloc(dataMatrix->size1, nFunctions);
    gsl_matrix *R = gsl_matrix_alloc(dataMatrix->size1, nFunctions);
    gsl_vector_view y = gsl_matrix_column(dataMatrix, 1);
    gsl_matrix *covCopy = gsl_matrix_alloc(covMatrix->size1, covMatrix->size2);
    gsl_matrix *covR = gsl_matrix_alloc(covMatrix->size1, covMatrix->size2);

    for (int i = 0; i < (A->size1); i++)
    {
        double x_i = gsl_matrix_get(dataMatrix, i, 0); //Gets the i'th x-value and stores it in x_i
        for (int j = 0; j < (A->size2); j++)
        {
            double A_ij = functions[j](x_i); //Computes f_j(x_i) and stores it in A_ij
            gsl_matrix_set(A, i, j, A_ij); //Set (i,j) in A equal to A_ij
        }
    }

    GS_decomp(A, R); //Gram-schmidt storing orthogonal Q in A and R is upper triangular
    GS_solve(A, R, (&y.vector), c); //Solves R*c = Q^t*b and stores in c

    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., R, R, 0, covCopy); //Compute R^t*R and store result in covCopy
    GS_decomp(covCopy, covR);
    GS_inverse(covCopy, covR, covMatrix); //Store inverse in covMatrix

    gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_matrix_free(covCopy);
    gsl_matrix_free(covR);
}

int main()
{
    int nData = 9;
    int nFunctions = 2;
    double x[] = {1., 2., 3., 4., 6., 9., 10., 13., 15.};
    double y[] = {117., 100., 88., 72., 53., 29.5, 25.2, 15.2, 11.1};

    gsl_matrix *dataMatrix = gsl_matrix_alloc(nData, 3);
    gsl_matrix *covariance = gsl_matrix_alloc(nFunctions, nFunctions);
    gsl_vector *c = gsl_vector_alloc(nFunctions);
    for (int i = 0; i < nData; i++)
    {
        gsl_matrix_set(dataMatrix, i, 0, x[i]); //Set first column in dataMatrix equal to x-values
        gsl_matrix_set(dataMatrix, i, 1, log(y[i])); //Set second column in dataMatrix equal to log of y-values
        gsl_matrix_set(dataMatrix, i, 2, (y[i] / 20) / y[i]); //Set third column in dataMatrix equal to uncertainty in y
    }

    f functions[2] = {&f_0, &f_1};


    least_sqr_fit(dataMatrix, functions, nFunctions, c, covariance);

    double c_0 = gsl_vector_get(c, 0);
    double c_1 = gsl_vector_get(c, 1);

    double dc_0 = sqrt(gsl_matrix_get(covariance, 0, 0));
    double dc_1 = sqrt(gsl_matrix_get(covariance, 1, 1));

    FILE *fitValues = fopen("out.fitvaluesexerciseAandB.txt", "w");

    fprintf(fitValues, "Parameters and half-life:\n");
    fprintf(fitValues, "c_0 = %g \n c_1 lambda = %g \n Calculated T_1/2 = %g +/- %g \n Real T_1/2 = %g \n", c_0, c_1,
            log(2) / (-c_1), dc_1, 3.63);
    fprintf(fitValues, "Difference between calculated half-life and modern value \n %g - %g = %g \n", log(2) / (-c_1),
            3.63, log(2) / (-c_1) - 3.63);
    fprintf(fitValues, "This value is not within the estimated uncertainty \n");
    fprintf(fitValues, "Covariance Matrix: \n");
    matrix_print(covariance, fitValues);

    //Fill txt file with data for datapoints
    FILE *data = fopen("out.dataplot.txt", "w");

    for (int i = 0; i < dataMatrix->size1; i++)
    {
        double x_i;
        double y_i;
        double dy_i;
        x_i = gsl_matrix_get(dataMatrix, i, 0);
        y_i = gsl_matrix_get(dataMatrix, i, 1);
        dy_i = gsl_matrix_get(dataMatrix, i, 2);
        fprintf(data, "%g %g %g\n", x_i, y_i, dy_i);
    }

    //Fill txt file with data for plot
    FILE *fit = fopen("out.fitplot.txt", "w");

    int resolutionForPlot = 200; // amount of points for plot

    for (int i = 0; i < resolutionForPlot; i++)
    {
        double x_i = ((double) i) / resolutionForPlot * x[nData - 1];
        double f_i = f_0(x_i) * c_0 + f_1(x_i) * c_1;
        double f_i_upper_bound = f_0(x_i) * (c_0 + dc_0) + f_1(x_i) * (c_1 + dc_1);
        double f_i_lower_bound = f_0(x_i) * (c_0 - dc_0) + f_1(x_i) * (c_1 - dc_1);
        fprintf(fit, "%g %g %g %g \n", x_i, f_i, f_i_upper_bound, f_i_lower_bound);
    }

    printf("\nResults of exercise A and B can be seen in file out.fitvaluesexerciseAandB.txt\n");
    printf("\nPlot from exercise A and C can be seen in fitplot.png\n\n");
    return 0;
}