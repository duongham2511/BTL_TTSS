#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "MatOp.h"

int main(int argc, char **argv)
{
    // int size_A = 1024;
    // int** matA = AllocateMemory2D<int>(size_A,size_A);
    // for (int i = 0; i <size_A; i ++)
    // {
    //     for (int j = 0; j < size_A; j++)
    //     {
    //         matA[i][j] = i+j;
    //     }
    // }
    // int size_B = 1024;
    // int** matB = AllocateMemory2D<int>(size_B,size_B);
    // for (int i = 0; i <size_B; i++)
    // {
    //     for (int j = 0; j <size_B; j++)
    //     {
    //         matB[i][j]= size_B-i-j;
    //     }
    // }

    // printMat(matA,size_A,size_A);
    // printf("\n");
    // printMat(matB,size_B,size_B);
    // printf("\n");

    // int size_C = 1024;
    // int** matC = AllocateMemory2D<int>(size_C,size_C);

    // double start_time_naive = omp_get_wtime();
    // matMul_Naive(matA,matB,matC,size_A);
    // double end_time_naive = omp_get_wtime();

    // double start_time_strassen = omp_get_wtime();
    // matMul_Strassen<int>(matA,matB,matC,size_A);
    // double end_time_strassen = omp_get_wtime();
    // // printMat(matC,size_C,size_C);
    // printf("Naive runtime: %f seconds\n",end_time_naive - start_time_naive);
    // printf("Strassen runtime: %f seconds\n",end_time_strassen - start_time_strassen);

    // FreeMemory2D<int>(matA);
    // FreeMemory2D<int>(matB);
    // FreeMemory2D<int>(matC);

    int **matA = AllocateMemory2D<int>(5,5);
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            matA[i][j] = i+j;
    // printMat(matA, 5,5);
    int new_size = pad2D<int>(matA,5);
    printf("%d\n",new_size);
    printMat(matA,new_size,new_size);

    FreeMemory2D<int>(matA);
    return 0;
}