#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_linalg.h>

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


double diffClock(clock_t startTime, clock_t endTime)
{
    /*
    Function diffClock, takes two arguments, the time stamp
    of the start and end of a time interval of interest.

    ¤ startTime : Timestamp of beginning of time interval
    ¤ endTime   : Timestamp of ending of time interval

    */

    double diffTicks = startTime - endTime;
    double diffms = (diffTicks * 10) / CLOCKS_PER_SEC;

    return diffms;
}


void timesJ(gsl_matrix *A, int p, int q, double theta)
{

    //Multiplies the matrix A with the Jacobi matrix J(p,q,theta) from the right

    double c = cos(theta);
    double s = sin(theta);

    for (int i = 0; i < (A->size1); i++)
    {
        double AJ_ip = c * gsl_matrix_get(A, i, p) - s * gsl_matrix_get(A, i, q);
        double AJ_iq = s * gsl_matrix_get(A, i, p) + c * gsl_matrix_get(A, i, q);

        gsl_matrix_set(A, i, p, AJ_ip);
        gsl_matrix_set(A, i, q, AJ_iq);
    }
}

void Jtimes(gsl_matrix *A, int p, int q, double theta)
{

    //Multiplies the matrix A with the Jacobi matrix J(p,q,theta) from the left

    double c = cos(theta);
    double s = sin(theta);

    for (int j = 0; j < (A->size2); j++)
    {
        double JA_pj = c * gsl_matrix_get(A, p, j) + s * gsl_matrix_get(A, q, j);
        double JA_qj = -s * gsl_matrix_get(A, p, j) + c * gsl_matrix_get(A, q, j);

        gsl_matrix_set(A, p, j, JA_pj);
        gsl_matrix_set(A, q, j, JA_qj);
    }
}


void jacobi_diag(gsl_matrix *A, gsl_matrix *V)
{

    // Jacobi eigenvalue algorithm for real symmetric matrices using cyclic sweeps. Convergence criterion is no change in eigenvalues after sweep.
    // A is turned into D, diagonal matrix with eigenvalues

    int changed = 0;
    int count = 0;
    int n = (A->size1);

    do
    {
        changed = 0;
        count++;
        for (int p = 0; p < n; p++)
        {
            for (int q = 0; q < n; q++)
            {
                if (p != q)
                {
                    double A_pq = gsl_matrix_get(A, p, q);
                    double A_pp = gsl_matrix_get(A, p, p);
                    double A_qq = gsl_matrix_get(A, q, q);
                    double theta = 0.5 * atan(2 * A_pq / (A_qq - A_pp));
                    double c = cos(theta);
                    double s = sin(theta);
                    double new_A_pp = c * c * A_pp - 2 * s * c * A_pq + s * s * A_qq;
                    double new_A_qq = s * s * A_pp + 2 * s * c * A_pq + c * c * A_qq;
                    if (new_A_pp != A_pp || new_A_qq != A_qq)
                    { // do rotation
                        changed = 1;
                        timesJ(A, p, q, theta);
                        Jtimes(A, p, q, -theta); // A←J^T*A*J
                        timesJ(V, p, q, theta); // V←V*J
                    }
                }
                else
                {}
            }
        }
    }
    while (changed != 0);
}

void jacobi_diag_ut(gsl_matrix *A, gsl_matrix *V)
{

    // Jacobi eigenvalue algorithm for real symmetric matrices using cyclic sweeps. Convergence criterion is no change in eigenvalues after sweep.
    // A is turned into D, diagonal matrix with eigenvalues. Faster runtime by considering only upper half of the matrix

    int changed = 0;
    int count = 0;
    int n = (A->size1);

    do
    {
        changed = 0;
        count++;
        for (int p = 0; p < n - 1; p++)
        {
            for (int q = p + 1; q < n; q++)
            {
                double A_pq = gsl_matrix_get(A, p, q);
                double A_pp = gsl_matrix_get(A, p, p);
                double A_qq = gsl_matrix_get(A, q, q);
                double theta = 0.5 * atan(2 * A_pq / (A_qq - A_pp));
                double c = cos(theta);
                double s = sin(theta);
                double new_A_pp = c * c * A_pp - 2 * s * c * A_pq + s * s * A_qq;
                double new_A_qq = s * s * A_pp + 2 * s * c * A_pq + c * c * A_qq;
                if (new_A_pp != A_pp || new_A_qq != A_qq)
                { // do rotation
                    changed = 1;
                    timesJ(A, p, q, theta);
                    Jtimes(A, p, q, -theta); // A←J^T*A*J
                    timesJ(V, p, q, theta); // V←V*J
                }
            }
        }
    }
    while (changed != 0);
}

int main()
{

    // Exercise A
    int nSize = 10; //Size of matrix
    FILE *exerciseA = fopen("out.exerciseA.txt", "w");

    //Allocating memory for matrices
    gsl_matrix *A = gsl_matrix_alloc(nSize, nSize);
    gsl_matrix *V = gsl_matrix_alloc(nSize, nSize);
    gsl_matrix *ACopy = gsl_matrix_alloc(nSize, nSize);
    gsl_matrix *VCopy = gsl_matrix_alloc(nSize, nSize);
    gsl_matrix *ATimesV = gsl_matrix_alloc(nSize, nSize);
    gsl_matrix *D = gsl_matrix_alloc(nSize, nSize);
    gsl_matrix *DCopy = gsl_matrix_alloc(nSize, nSize);
    gsl_matrix *VtimesV = gsl_matrix_alloc(nSize, nSize);
    gsl_matrix *DtimesV = gsl_matrix_alloc(nSize, nSize);

    //Generating random symmetric matrix
    for (int i = 0; i < nSize; i++)
    {
        for (int j = i; j < nSize; j++)
        {
            double A_ij = ((double) rand() / RAND_MAX);
            gsl_matrix_set(A, i, j, A_ij);
            gsl_matrix_set(A, j, i, A_ij);
        }
    }

    gsl_matrix_set_identity(V); //Turn V into identity matrix
    gsl_matrix_memcpy(ACopy, A); //Turn ACopy into copy of A

    fprintf(exerciseA, "Symmetric matrix A = \n");
    matrix_print(A, exerciseA);

    jacobi_diag(A, V); //Do Jacobi eigenvalue algorithm. Turns A into D

    gsl_matrix_memcpy(D, A); //Turn D into copy of A (Which was turned into D by jacobi diagonalization)
    gsl_matrix_memcpy(DCopy, D); //Turn DCopy into copy of D

    fprintf(exerciseA, "\n Diagonal matrix with eigenvalues D = \n");
    matrix_print(A, exerciseA);

    fprintf(exerciseA, "\n Is V^T*A*V = D? We check if V^T*A*V-D = 0 \n");
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ACopy, V, 0., ATimesV); //Computes A*V and stores in ATimesV
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, ATimesV, -1., D); //Computes V^T*A*V-D and stores in D
    matrix_print(D, exerciseA);

    fprintf(exerciseA, "\n Is V*D*V^T = A? We check if V*D*V^T-A =0 \n");
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DCopy, V, 0., DtimesV); //Computes D*V^T and stores in DtimesV
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, DtimesV, -1., ACopy); //Computes V*D*V^T - D and stores in Dcopy
    matrix_print(ACopy, exerciseA);

    fprintf(exerciseA, "\n Is V^T*V = 1? \n");
    gsl_matrix_memcpy(VCopy, V); //Turns Vcopy into copy of V
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, VCopy, V, 0., VtimesV);
    matrix_print(VtimesV, exerciseA);

    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_matrix_free(ACopy);
    gsl_matrix_free(VCopy);
    gsl_matrix_free(ATimesV);
    gsl_matrix_free(D);
    gsl_matrix_free(DCopy);
    gsl_matrix_free(VtimesV);
    gsl_matrix_free(DtimesV);

    //Exercise B

    FILE *exerciseB = fopen("out.exerciseB.txt", "w");

    int hamilSize = 20; //Size of hamilton matrix
    double s = 1 / (double) (hamilSize + 1);
    gsl_matrix *H = gsl_matrix_alloc(hamilSize, hamilSize);
    for (int i = 0; i < hamilSize - 1; i++)
    {
        gsl_matrix_set(H, i, i, -2);
        gsl_matrix_set(H, i, i + 1, 1);
        gsl_matrix_set(H, i + 1, i, 1);
    }
    gsl_matrix_set(H, hamilSize - 1, hamilSize - 1, -2);
    gsl_matrix_scale(H, -1 / s / s);

    gsl_matrix *W = gsl_matrix_alloc(hamilSize, hamilSize);
    gsl_matrix_set_identity(W);

    jacobi_diag(H, W); //Turns H into diagonal matrix with energies

    //Dmitri's method seems to not be sorted correctly? We sort the eigenvalues in the diagonal
    gsl_vector_view energiesH = gsl_matrix_diagonal(H);
    gsl_vector *indexForSort = gsl_vector_alloc(hamilSize);
    gsl_vector *energiesHCopy = gsl_vector_alloc(hamilSize);
    gsl_vector_memcpy(energiesHCopy, &energiesH.vector);

    for (int i = 0; i < hamilSize; i++)
    {
        gsl_vector_set(indexForSort, i, i);
    }

    gsl_sort_vector2(energiesHCopy, indexForSort);

    fprintf(exerciseB, "Energies for quantum particle in a box. Second column is calculated, third is exact: \n");
    for (int k = 0; k < hamilSize / 3; k++)
    {
        double exact = M_PI * M_PI * (k + 1) * (k + 1);
        double calculated = gsl_vector_get(energiesHCopy, k);
        fprintf(exerciseB, "%i %g %g\n", k, calculated, exact);
    }

    FILE *psiPlot = fopen("out.psiPlot.txt", "w");

    fprintf(psiPlot, "%g %g %g %g %g %g %g\n", 0., 0., 0., 0., 0., 0., 0.);
    for (int k = 0; k < hamilSize; k++)
    {
        int a_1 = gsl_vector_get(indexForSort,
                                 0); //We need to pull out the eigenvectors corresponding to the eigenvalues (lowest energy) So we keep track of the index vector to find right eigenvectors
        int a_2 = gsl_vector_get(indexForSort, 1);
        int a_3 = gsl_vector_get(indexForSort, 2);
        fprintf(psiPlot, "%g %g %g %g %g %g %g\n", (k + 1) / (double) (hamilSize + 1), 3.2 * gsl_matrix_get(W, k, a_1),
                3.2 * gsl_matrix_get(W, k, a_2), 3.2 * gsl_matrix_get(W, k, a_3), sin(s * M_PI * (k + 1)),
                sin(2 * s * M_PI * (k + 1)), -sin(3 * s * M_PI * (k + 1)));
    } //The scaling used above has been found by simply looking at the plot. Since the purpose is to compare, the exact size does not matter.
    fprintf(psiPlot, "%g %g %g %g %g %g %g\n", 1., 0., 0., 0., 0., 0., 0.);

    gsl_matrix_free(H);
    gsl_matrix_free(W);
    gsl_vector_free(indexForSort);
    gsl_vector_free(energiesHCopy);

    //Exercise C


    int firstLoop = 0;
    double baseTime = 0;

    FILE *jacobi_timer = fopen("out.jacobi_timer.txt", "w");
    FILE *testfile = fopen("out.test.txt", "w");


    for (int q = 40; q < 80; q += 2)
    {

        gsl_matrix *randomA = gsl_matrix_alloc(q, q);
        gsl_matrix *randomV = gsl_matrix_alloc(q, q);
        gsl_matrix *randomUpperA = gsl_matrix_alloc(q, q);
        gsl_matrix *randomUpperV = gsl_matrix_alloc(q, q);
        gsl_matrix *randomGslA = gsl_matrix_alloc(q, q);
        gsl_matrix *randomGslV = gsl_matrix_alloc(q, q);
        gsl_vector *gslS = gsl_vector_alloc(q);


        for (int i = 0; i < q; i++)
        { //Create random symmetric matrix
            for (int j = i; j < q; j++)
            {
                double A_ij = ((double) rand() / RAND_MAX);
                gsl_matrix_set(randomA, i, j, A_ij);
                gsl_matrix_set(randomA, j, i, A_ij);
            }
        }
        gsl_matrix_set_identity(randomV);

        gsl_matrix_memcpy(randomUpperA, randomA);
        gsl_matrix_memcpy(randomUpperV, randomV);
        gsl_matrix_memcpy(randomGslA, randomA);
        gsl_matrix_memcpy(randomGslV, randomV);


        clock_t beginMine = clock();
        clock_t endMine = clock();
        clock_t beginUT = clock();
        clock_t endUT = clock();
        clock_t beginGSL = clock();
        clock_t endGSL = clock();

        matrix_print(randomA, testfile);

        beginMine = clock();
        jacobi_diag(randomA, randomV);
        endMine = clock();

        beginUT = clock();
        jacobi_diag_ut(randomUpperA, randomUpperV);
        endUT = clock();

        beginGSL = clock();
        gsl_linalg_SV_decomp_jacobi(randomGslA, randomGslV, gslS);
        endGSL = clock();

        if (firstLoop == 0)
        {
            baseTime += (double) (diffClock(endMine, beginMine));
            firstLoop = 1;
        }

        double O3 = pow(((double) q) / (40), 3) * baseTime;
        fprintf(jacobi_timer, "%d\t%g\t%g\t%g\t%g\n", q, (double) (diffClock(endMine, beginMine)), O3,
                (double) (diffClock(endUT, beginUT)), (double) (diffClock(endGSL, beginGSL)));


        gsl_matrix_free(randomA);
        gsl_matrix_free(randomV);
        gsl_matrix_free(randomUpperA);
        gsl_matrix_free(randomUpperV);
        gsl_matrix_free(randomGslA);
        gsl_matrix_free(randomGslV);
        gsl_vector_free(gslS);

        printf("Done with %d dimension matrix\n", q);


    }

    printf("\n Result of exercise A can be seen in file out.exerciseA.txt \n\n");
    printf("Result of exercise B can be seen in file out.exerciseB.txt \n plot of wavefunctions in exercise B can be seen in psiplot.png \n\n");
    printf("Result of exercise C can be seen in compareplot.png\n\n");

    return 0;
}