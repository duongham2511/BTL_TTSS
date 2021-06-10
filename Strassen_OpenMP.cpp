#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "MatOp.h"

#define SIZE 4096

int main(int argc, char **argv)
{
    const int n = SIZE;
    double **A = AllocateMemory2D<double>(n,n);
    double **B = AllocateMemory2D<double>(n,n);
    double **C = AllocateMemory2D<double>(n,n);

    int i,j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = (double)i / (double)n;
            B[i][j] = (double)j / (double)n;
        }
    }

    matMul_Naive(A,B,C,n);

    double norm = 0.0;
    for (i=0  ; i < n ; i++)
        for (j=0  ; j < n ; j++)
            norm += (C[i][j]-(double)(i*j)/(double)n)*(C[i][j]-(double)(i*j)/(double)n);

    if (norm > 1e-10)
        printf("Something is wrong... Norm is equal to %f\n", norm);
    else
        printf("Yup, we're good!\n");

    FreeMemory2D<double>(A);
    FreeMemory2D<double>(B);
    FreeMemory2D<double>(C);
    return 0;
}