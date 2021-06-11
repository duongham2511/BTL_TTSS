#ifndef TESTSUITE_H
#define TESTSUITE_H

#include <omp.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "Utils.h"
#include "MatOp.h"
#include "MatOp_P.h"
#include "Strassen.h"

class TestSuite{

    static void GenTestMat(double** A, double** B, int rowA, int colA_rowB, int colB)
    {
        int i,j;
        for (j = 0; j < colA_rowB; j++) {
            for (i = 0; i < rowA; i++)
            {
                A[i][j] = (double)i/double(colA_rowB);
            }
            for (i = 0; i < colB; i++)
            {
                B[j][i] = (double)i/double(colA_rowB);
            }
        }
    }

    static bool TestProduct(double** C, int rowA, int colA_rowB, int colB)
    {
        int i,j;
        double norm = 0.0;
        for (i=0  ; i < rowA ; i++)
            for (j=0  ; j < colB ; j++)
                norm += (C[i][j]-(double)(i*j)/(double)colA_rowB)*(C[i][j]-(double)(i*j)/(double)colA_rowB);

        if (norm > 1e-10)
        {
            printf("Something is wrong... Norm is equal to %f\n", norm);
            return false;
        } else
        {
            printf("Yup, we're good!\n");
            return true;
        }
    }

    public:
    static void NaiveSerialTest(int rowA, int colA_rowB, int colB)
    {
        double **A = Utility::AllocateMemory2D<double>(rowA,colA_rowB);
        double **B = Utility::AllocateMemory2D<double>(colA_rowB,colB);
        double **C = Utility::AllocateMemory2D<double>(rowA,colB);
        TestSuite::GenTestMat(A,B,rowA,colA_rowB,colB);

        printf("Multiplying a %d x %d matrix with a %d x %d matrix using naive serial.\n",rowA, colA_rowB, colA_rowB, colB);
        double start = omp_get_wtime();
        MatOp::matMul_Naive(A,B,C,rowA,colA_rowB,colB);
        double end = omp_get_wtime();

        if (TestSuite::TestProduct(C,rowA,colA_rowB,colB))
            printf("Runtime is %.2f seconds.\n",end - start);

        Utility::FreeMemory2D<double>(A);
        Utility::FreeMemory2D<double>(B);
        Utility::FreeMemory2D<double>(C);
    }

    static void NaiveOMPTest(int rowA, int colA_rowB, int colB)
    {
        double **A = Utility::AllocateMemory2D<double>(rowA,colA_rowB);
        double **B = Utility::AllocateMemory2D<double>(colA_rowB,colB);
        double **C = Utility::AllocateMemory2D<double>(rowA,colB);
        TestSuite::GenTestMat(A,B,rowA,colA_rowB,colB);

        printf("Multiplying a %d x %d matrix with a %d x %d matrix using naive OMP.\n",rowA, colA_rowB, colA_rowB, colB);
        double start = omp_get_wtime();
        MatOp_P::matMul_Naive(A,B,C,rowA,colA_rowB,colB);
        double end = omp_get_wtime();

        if (TestSuite::TestProduct(C,rowA,colA_rowB,colB))
            printf("Runtime is %.2f seconds.\n",end - start);

        Utility::FreeMemory2D<double>(A);
        Utility::FreeMemory2D<double>(B);
        Utility::FreeMemory2D<double>(C);
    }

    static void StrassenSerialTest(int rowA, int colA_rowB, int colB)
    {
        double **A = Utility::AllocateMemory2D<double>(rowA,colA_rowB);
        double **B = Utility::AllocateMemory2D<double>(colA_rowB,colB);
        double **C;
        TestSuite::GenTestMat(A,B,rowA,colA_rowB,colB);

        printf("Multiplying a %d x %d matrix with a %d x %d matrix using Strassen Serial.\n",rowA, colA_rowB, colA_rowB, colB);

        int max = fmax(rowA,colB);
        max = max > colA_rowB ? max : colA_rowB;

        int n = Utility::getPower2(max);
        bool needPad = !((n == max) && (rowA == colA_rowB) && (colA_rowB == colB));
        if  (needPad)
        {
            printf("Matrix A of %d x %d padded to %d x %d.\n", rowA, colA_rowB, n, n);
            MatOp::pad2D(A,rowA,colA_rowB,n,n);
            // Utility::printMat(A,n,n);
            printf("Matrix B of %d x %d padded to %d x %d.\n", colA_rowB, colB, n, n);
            MatOp::pad2D(B,colA_rowB,colB,n,n);
            // Utility::printMat(B,n,n);
            needPad = true;
        }

        C = Utility::AllocateMemory2D<double>(n,n);
        double start = omp_get_wtime();
        matMul_Strassen_S(A,B,C,n);
        double end = omp_get_wtime();

        if (needPad)
        {
            printf("Matrix C of %d x %d shrunk to %d x %d.\n", n, n, rowA, colB);
            MatOp::shrink2D<double>(C,n,n,rowA,colB);
            // Utility::printMat(C,rowA,colB);
        }

        if (TestSuite::TestProduct(C,rowA,colA_rowB,colB))
            printf("Runtime is %.2f seconds.\n",end - start);

        Utility::FreeMemory2D<double>(A);
        Utility::FreeMemory2D<double>(B);
        Utility::FreeMemory2D<double>(C);
    }
};

#endif