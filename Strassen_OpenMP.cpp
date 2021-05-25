#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "MatOp.h"

#define SIZE 2

int main(int argc, char **argv)
{
    int size_A = SIZE;
    int** matA = AllocateMemory2D<int>(size_A,size_A);
    for (int i = 0; i <size_A; i ++)
    {
        for (int j = 0; j < size_A; j++)
        {
            matA[i][j] = i+j;
        }
    }
    int size_B = SIZE;
    int** matB = AllocateMemory2D<int>(size_B,size_B);
    for (int i = 0; i <size_B; i++)
    {
        for (int j = 0; j <size_B; j++)
        {
            matB[i][j]= size_B-i-j;
        }
    }

    printMat(matA,size_A,size_A);
    printf("\n");
    printMat(matB,size_B,size_B);
    printf("\n");

    int size_C = SIZE;
    int** matC = AllocateMemory2D<int>(size_C,size_C);


    // double start_time_naive = omp_get_wtime();
    // matMul_Naive(matA,matB,matC,size_A);
    // double end_time_naive = omp_get_wtime();

    // double start_time_strassen = omp_get_wtime();
    matMul_Strassen_v2<int>(matA,matB,matC,size_A);
    // double end_time_strassen = omp_get_wtime();
    printMat(matC,size_C,size_C);
    // printf("Naive runtime: %f seconds\n",end_time_naive - start_time_naive);
    // printf("Strassen runtime: %f seconds\n",end_time_strassen - start_time_strassen);

    FreeMemory2D<int>(matA);
    FreeMemory2D<int>(matB);
    FreeMemory2D<int>(matC);

    // int **matA = AllocateMemory2D<int>(8,8);
    // for (int i = 0; i < 8; i++)
    //     for (int j = 0; j < 8; j++)
    //         matA[i][j] = i+j;
    // // printMat(matA, 8,8);
    // int **matA11 = AllocateMemory2D<int>(4,4);
    // int **matA22 = AllocateMemory2D<int>(4,4);
    // split2D<int>(matA,matA22,4,4,4);
    // split2D<int>(matA,matA11,4,0,0);
    // join2D<int>(matA11,matA,4,4,4);
    // join2D<int>(matA22,matA,4,0,0);
    // printMat(matA,8,8);

    // FreeMemory2D<int>(matA);
    // FreeMemory2D<int>(matA11);
    // FreeMemory2D<int>(matA22);
    return 0;
}