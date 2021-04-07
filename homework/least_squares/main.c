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

void least_sqr_fit(gsl_matrix* dataMatrix, f functions[],int nFunctions, gsl_vector* c, gsl_matrix* covMatrix){
    gsl_matrix* A = gsl_matrix_alloc(dataMatrix->size1,nFunctions);
    gsl_matrix* R = gsl_matrix_alloc(dataMatrix->size1,nFunctions);
    gsl_vector_view y = gsl_matrix_column(dataMatrix,1);
    gsl_matrix* covCopy = gsl_matrix_alloc(covMatrix->size1, covMatrix->size2);
    gsl_matrix* covR = gsl_matrix_alloc(covMatrix->size1,covMatrix->size2);

    for(int i=0; i<(A->size1); i++){
        double x_i = gsl_matrix_get(dataMatrix,i,0); //Gets the i'th x-value and stores it in x_i
        for (int j=0; j<(A->size2); j++){
            double A_ij = functions[j](x_i); //Computes f_j(x_i) and stores it in A_ij
            gsl_matrix_set(A,i,j,A_ij); //Set (i,j) in A equal to A_ij
        }
    }

    GS_decomp(A,R); //Gram-schmidt storing orthogonal Q in A and R is upper triangular
    GS_solve(A,R,(&y.vector),c); //Solves R*c = Q^t*b and stores in c

    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,R,R,0,covCopy); //Compute R^t*R and store result in covCopy
    GS_decomp(covCopy,covR);
    GS_inverse(covCopy,covR,covMatrix); //Store inverse in covMatrix

    gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_matrix_free(covCopy);
    gsl_matrix_free(covR);
}

int main(){
    int nData = 9;
    int nFunctions = 2;
    double x[] = {1., 2., 3., 4., 6., 9., 10., 13., 15.};
    double y[] = {117., 100., 88., 72.,  53., 29.5, 25.2, 15.2, 11.1};

    gsl_matrix* dataMatrix = gsl_matrix_alloc(nData,3);
    gsl_matrix* covariance = gsl_matrix_alloc(nFunctions,nFunctions);
    gsl_vector* c = gsl_vector_alloc(nFunctions);
    for(int i=0; i<nData; i++){
        gsl_matrix_set(dataMatrix,i,0,x[i]); //Set first column in dataMatrix equal to x-values
        gsl_matrix_set(dataMatrix,i,1,log(y[i])); //Set second column in dataMatrix equal to log of y-values
        gsl_matrix_set(dataMatrix,i,2,(y[i]/20)/y[i]); //Set third column in dataMatrix equal to uncertainty in y
    }

    f functions[2] = {&f_0, &f_1};


    least_sqr_fit(dataMatrix,functions,nFunctions,c,covariance);

    double c_0 = gsl_vector_get(c,0);
    double c_1 = gsl_vector_get(c,1);

    double dc_0 = sqrt(gsl_matrix_get(covariance,0,0));
    double dc_1 = sqrt(gsl_matrix_get(covariance,1,1));

    FILE* fitValues = fopen("out.fitvaluesexerciseAandB.txt","w");

    fprintf(fitValues,"Parameters and half-life:\n");
    fprintf(fitValues,"c_0 = %g \n c_1 lambda = %g \n Calculated T_1/2 = %g +/- %g \n Real T_1/2 = %g \n", c_0,c_1,log(2)/(-c_1),dc_1,3.63);
    fprintf(fitValues, "Difference between calculated half-life and modern value \n %g - %g = %g \n",log(2)/(-c_1),3.63,log(2)/(-c_1)-3.63);
    fprintf(fitValues, "This value is not within the estimated uncertainty \n");
    fprintf(fitValues,"Covariance Matrix: \n");
    matrix_print(covariance,fitValues);

    //Fill txt file with data for datapoints
    FILE* data = fopen("out.dataplot.txt","w");

    for(int i=0; i<dataMatrix->size1; i++){
        double x_i;
        double y_i;
        double dy_i;
        x_i=gsl_matrix_get(dataMatrix,i,0);
        y_i=gsl_matrix_get(dataMatrix,i,1);
        dy_i=gsl_matrix_get(dataMatrix,i,2);
        fprintf(data,"%g %g %g\n",x_i,y_i,dy_i);
    }

    //Fill txt file with data for plot
    FILE* fit = fopen("out.fitplot.txt","w");

    int resolutionForPlot = 200; // amount of points for plot

    for(int i=0; i<resolutionForPlot; i++){
        double x_i = ((double) i)/resolutionForPlot*x[nData-1];
        double f_i = f_0(x_i)*c_0+f_1(x_i)*c_1;
        double f_i_upper_bound = f_0(x_i)*(c_0+dc_0)+f_1(x_i)*(c_1+dc_1);
        double f_i_lower_bound = f_0(x_i)*(c_0-dc_0)+f_1(x_i)*(c_1-dc_1);
        fprintf(fit,"%g %g %g %g \n",x_i,f_i,f_i_upper_bound,f_i_lower_bound);
    }

    return 0;
}