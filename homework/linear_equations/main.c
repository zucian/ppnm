#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <time.h>

double randomNumber( unsigned int *seed ){
    double maxRand      =   (double)RAND_MAX;           // Maximum random number, cast to double
    double randNum      =   (double)rand_r( seed );     // Generate pseudo-random number from seed, cast to double
    return randNum/maxRand;
}

double diffClock(clock_t startTime, clock_t endTime){
    /*
    Function diffClock, takes two arguments, the time stamp
    of the start and end of a time interval of interest.

    ¤ startTime : Timestamp of beginning of time interval
    ¤ endTime   : Timestamp of ending of time interval

    */

    double diffTicks  =  startTime - endTime;
    double diffms     =  (diffTicks * 10) / CLOCKS_PER_SEC;

    return diffms;
}

void set_data_tall(gsl_matrix* tallMatrix, gsl_vector* tallRHSVector, unsigned int *seed){
    /* Method to set the data of a tall matrix using pseudorandom numbers.
     *
     *  ¤ gsl_matrix*    tallMatrix     : An arbitrary n x m matrix (n >= m)
     *  ¤ gsl_vector*    tallRHSVector  : A vector of dimension n x 1
     *  ¤ unsigned int*  seed           : A seed for the pseudorandom number generator
     *
     */

    for(int rowId = 0; rowId < tallMatrix -> size1; rowId++){
        for(int colId = 0; colId < tallMatrix -> size2; colId++){
            gsl_matrix_set(tallMatrix, rowId, colId, randomNumber(seed));
        }
        gsl_vector_set(tallRHSVector, rowId, randomNumber(seed));
    }
}

void vector_print(gsl_vector* v, FILE* file){
    for(int i=0;i<v->size;i++)fprintf(file, "%10g ",gsl_vector_get(v,i));
}

/*
void matrix_print(int numOfRows, gsl_matrix* matrixToPrint, char* string ){
    printf("\n%s\n", string);
    for (int rowId = 0; rowId < numOfRows; rowId++){
        gsl_vector_view matrixToPrint_row = gsl_matrix_row (matrixToPrint, rowId);
        gsl_vector* vector = &matrixToPrint_row.vector;
        for(int iter = 0; iter < vector -> size; iter++){
            if ( gsl_vector_get(vector, iter) > 1e-10 ){
                printf("%10g\t", gsl_vector_get(vector, iter));
            }
            else { printf("%10g\t", 0.0); }
        }
        printf("\n");
    }
}
 */

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
    gsl_vector* a_i = gsl_vector_alloc(A->size1); //Allocate memory for a_i and q_i to be used later.
    gsl_vector* q_i = gsl_vector_alloc(A->size1); //Size should be equal to amount of rows in A

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

    //Takes a upper triangular matrix R and performs back substitution. Returns solution in x.

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

    gsl_blas_dgemv(CblasTrans, 1., Q,b,0.,x); //Computes the matrix vector product Q^T*b and saves the result in x
    backsub(R,x); //Compute back substitution of R and x and returns answer in x
}

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
    gsl_vector* unitVector = gsl_vector_alloc(R->size1);
    gsl_matrix* RInverse = gsl_matrix_alloc(R->size1,R->size2);
    for(int i=0; i<(R->size2); i++){ //Loops over every column in R
        gsl_vector_set_basis(unitVector,i);
        backsub(R,unitVector); //Solve R*x = e_i, where e_i is the i'th standard basis vector
        gsl_matrix_set_col(RInverse, i,unitVector); //Turns i'th column in R^-1 into solution to R*x = e_i
    }
    gsl_vector_free(unitVector);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,RInverse,Q,0.,B);
}

int main(){
    int n=10;
    int m=10;

    FILE* exerciseA = fopen("out.exerciseA.txt","w");


    gsl_matrix* A = gsl_matrix_alloc(n,m);
    gsl_matrix* A_copy = gsl_matrix_alloc(n,m);
    gsl_matrix* R = gsl_matrix_alloc(m,m);

    for(int i=0; i<(A->size1); i++){ //Creates a random matrix of size n,m
        for(int j=0; j<(A->size2); j++){
            double A_ij = rand()/100;
            gsl_matrix_set(A,i,j,A_ij);
        }
    }

    gsl_matrix_memcpy(A_copy,A); //Makes a copy of A in A_copy

    GS_decomp(A,R); //Performs gram-schmidt on A, stores Q in A and upper triangular in R

    gsl_matrix* Q_test = gsl_matrix_alloc(m,m); //Allocates memory for Q^T*Q matrix that tests orthogonality of Q
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,A,A,0.,Q_test); //Makes Q-test into Q^t*Q

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,A,R,-1.,A_copy); //Makes A_copy into Q*R-A

    fprintf(exerciseA,"Is R upper triangular? \n");
    matrix_print(R,exerciseA);

    fprintf(exerciseA,"\n Is Q^T*Q=1? \n");
    matrix_print(Q_test,exerciseA);

    fprintf(exerciseA,"\n Is QR = A? We compute Q*R-A and show it is equal to 0 \n");
    matrix_print(A_copy,exerciseA);

    gsl_vector* b = gsl_vector_alloc(n);
    gsl_vector* x = gsl_vector_alloc(m);

    for(int i=0; i<(A->size1); i++){ //Creates a random matrix of size n,m
        for(int j=0; j<(A->size2); j++){
            double A_ij = rand()/100;
            gsl_matrix_set(A,i,j,A_ij);
        }
    }

    for (int i=0; i<(b->size); i++){ //Creates a random vector of size n
        double b_i = rand()/100;
        gsl_vector_set(b,i,b_i);
    }

    gsl_matrix_memcpy(A_copy,A); //Makes copy of A and stores in A_copy

    GS_decomp(A,R); //Turns A into Q*R. Q is stored in A.

    GS_solve(A,R,b,x); //Solves Q*R*x = b and stores result in x

    gsl_blas_dgemv(CblasNoTrans,1.,A_copy,x,-1.,b); //Computes A*x-b and stores it in b

    //vector_print("Is A*x = b? We compute A*x-b and show it is equal to zero",b);

    fprintf(exerciseA, "\n Is A*x = b? We compute A*x-b and show it is equal to zero \n");
    vector_print(b,exerciseA);

    FILE* exerciseB = fopen("out.exerciseB.txt","w");

    for(int i=0; i<(A->size1); i++){ //Creates a random matrix of size n,m
        for(int j=0; j<(A->size2); j++){
            double A_ij = rand()/100;
            gsl_matrix_set(A,i,j,A_ij);
        }
    }

    gsl_matrix_memcpy(A_copy,A);
    gsl_matrix* B = gsl_matrix_alloc(n,m);
    gsl_matrix* I = gsl_matrix_alloc(n,m);

    GS_decomp(A,R);

    GS_inverse(A,R,B);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,A_copy,B,0.,I);

    fprintf(exerciseB,"Is A*B=I? \n");
    matrix_print(I,exerciseB);


    FILE* GS_timer = fopen("out.GS_timer.txt","w");
    /*
    int count = 0;
    int i = 100;

    while(count<100){
        gsl_matrix* A_time = gsl_matrix_alloc(i,i);
        gsl_matrix* R_time = gsl_matrix_alloc(i,i);
        gsl_matrix* gsl_A_time = gsl_matrix_alloc(i,i);
        gsl_vector* gsl_V_time = gsl_vector_alloc(i);

        for(int i=0; i<(A_time->size1); i++){ //Creates a random matrix of size n,m
            for(int j=0; j<(A_time->size2); j++){
                double A_ij = rand()/RAND_MAX;
                gsl_matrix_set(A_time,i,j,A_ij);
            }
        }

        gsl_matrix_memcpy(gsl_A_time,A_time);

        clock_t beginMine = clock();
        clock_t endMine   = clock();
        clock_t beginGSL = clock();
        clock_t endGSL = clock();

        beginMine = clock();
        GS_decomp(A_time,R_time);
        endMine = clock();

        double baseTime = 0;

        if(i=100){
            baseTime = (double)(diffClock(endMine,beginMine));
        }

        beginGSL = clock();
        gsl_linalg_QR_decomp(gsl_A_time,gsl_V_time);
        endGSL = clock();

        double O3 = pow(((double) i)/(100),3)*baseTime;
        fprintf(GS_timer, "%d\t%g\t%g\n", i, (double)(diffClock(endMine,beginMine)),O3);

        gsl_matrix_free(A_time);
        gsl_matrix_free(R_time);

        i += 1;
        count++;
    }
     */
    int firstLoop = 0;
    double baseTime = 0;

    for(int q=200; q<400; q++){
        gsl_matrix* A_time = gsl_matrix_alloc(q,q);
        gsl_matrix* R_time = gsl_matrix_alloc(q,q);
        gsl_matrix* gsl_A_time = gsl_matrix_alloc(q,q);
        gsl_vector* gsl_V_time = gsl_vector_alloc(q);

        for(int i=0; i<(A_time->size1); i++){ //Creates a random matrix of size n,m
            for(int j=0; j<(A_time->size2); j++){
                double A_ij = rand()/RAND_MAX;
                gsl_matrix_set(A_time,i,j,A_ij);
            }
        }

        gsl_matrix_memcpy(gsl_A_time,A_time);

        clock_t beginMine = clock();
        clock_t endMine   = clock();
        clock_t beginGSL = clock();
        clock_t endGSL = clock();

        beginMine = clock();
        GS_decomp(A_time,R_time);
        endMine = clock();

        if(firstLoop==0){
            baseTime += (double)(diffClock(endMine,beginMine));
            firstLoop = 1;
        }

        beginGSL = clock();
        gsl_linalg_QR_decomp(gsl_A_time,gsl_V_time);
        endGSL = clock();

        double O3 = pow(((double) q)/(200),3)*baseTime;
        fprintf(GS_timer, "%d\t%g\t%g\t%g\n", q, (double)(diffClock(endMine,beginMine)),O3,(double)(diffClock(endGSL,beginGSL)));

        gsl_matrix_free(A_time);
        gsl_matrix_free(R_time);
        gsl_matrix_free(gsl_A_time);
        gsl_vector_free(gsl_V_time);

        if(q%10==0){
            printf("Done with %d \n",q);
        }

    }


    return 0;
}