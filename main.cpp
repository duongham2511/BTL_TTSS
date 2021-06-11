#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "Utils.h"
#include "MatOp.h"
#include "MatOp_P.h"
#include "Strassen.h"
#include "TestSuite.h"

#define SIZE 4

int main(int argc, char **argv)
{
    // const int n = SIZE;
    omp_set_num_threads(8);
    TestSuite::StrassenSerialTest(300,500,200);
    // double **A = Utility::AllocateMemory2D<double>(n,n);
    // double **B = Utility::AllocateMemory2D<double>(n,n);
    // double **C = Utility::AllocateMemory2D<double>(n,n);

    // int i,j;

    // for (i = 0; i < n; i++) {
    //     for (j = 0; j < n; j++) {
    //         A[i][j] = (double)i / (double)n;
    //         B[i][j] = (double)j / (double)n;
    //     }
    // }

    // Utility::printMat(A,n,n);
    // printf("\n");
    // Utility::printMat(B,n,n);

    // double start_strassen_parallel = omp_get_wtime();
    // // matMul_Strassen_P(A,B,C,n);
    // MatOp::matMul_Naive(A,B,C,n,n,n);
    // printf("\n");
    // Utility::printMat(C,n,n);
    // double end_strassen_parallel = omp_get_wtime();

    // double norm = 0.0;
    // for (i=0  ; i < n ; i++)
    //     for (j=0  ; j < n ; j++)
    //         norm += (C[i][j]-(double)(i*j)/(double)n)*(C[i][j]-(double)(i*j)/(double)n);

    // if (norm > 1e-10)
    //     printf("Something is wrong... Norm is equal to %f\n", norm);
    // else
    // {
    //     printf("Yup, we're good!\n");
    //     printf("Runtime is %.2f seconds.\n",end_strassen_parallel - start_strassen_parallel);
    // }

    // Utility::FreeMemory2D<double>(A);
    // Utility::FreeMemory2D<double>(B);
    // Utility::FreeMemory2D<double>(C);
    return 0;
}