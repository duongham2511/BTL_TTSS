#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "MatOp.h"

int main(int argc, char **argv)
{
    int size_A = 4;
    int** matA = AllocateMemory2D<int>(size_A,size_A);
    for (int i = 0; i <size_A; i ++)
    {
        for (int j = 0; j < size_A; j++)
        {
            matA[i][j] = i+j;
        }
    }
    int size_B = 4;
    int** matB = AllocateMemory2D<int>(size_B,size_B);
    for (int i = 0; i <size_B; i++)
    {
        for (int j = 0; j <size_B; j++)
        {
            matB[i][j]= size_B-i-j;
        }
    }

    printMat(matA,size_A,size_A);
    printMat(matB,size_B,size_B);

    int size_C = 4;
    int** matC = AllocateMemory2D<int>(size_C,size_C);

    matMul_Strassen<int>(matA,matB,matC,size_A);
    printMat(matC,size_C,size_C);

    FreeMemory2D<int>(matA);
    FreeMemory2D<int>(matB);
    FreeMemory2D<int>(matC);

    return 0;
}