#ifndef STRASSEN_H
#define STRASSEN_H

#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "MatOp.h"
#include "Utils.h"

#define THRESHOLD 64

template <typename T>
void matMul_Strassen_S(T**& A, T**& B, T**& C, int size)
{
    if (size <= THRESHOLD)
    {
        MatOp::matMul_Naive<T>(A, B, C, size);
    } else 
    {
        T** T_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_4 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_5 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_6 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_7 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_8 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_9 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_10 = Utility::AllocateMemory2D<T>(size/2,size/2);

        T** B_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** B_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** B_12 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** B_21 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** A_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** A_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** A_12 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** A_21 = Utility::AllocateMemory2D<T>(size/2,size/2);

        T** Q_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_4 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_5 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_6 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_7 = Utility::AllocateMemory2D<T>(size/2,size/2);

        MatOp::split2D<T>(A,A_11,size/2,0,0);
        MatOp::split2D<T>(A,A_12,size/2,0,size/2);
        MatOp::split2D<T>(A,A_21,size/2,size/2,0);
        MatOp::split2D<T>(A,A_22,size/2,size/2,size/2);
        MatOp::split2D<T>(B,B_11,size/2,0,0);
        MatOp::split2D<T>(B,B_12,size/2,0,size/2);
        MatOp::split2D<T>(B,B_21,size/2,size/2,0);
        MatOp::split2D<T>(B,B_22,size/2,size/2,size/2);

        MatOp::MatOp::AddMat<T>(A_11,A_22,T_1,size/2);
        MatOp::AddMat<T>(B_11,B_22,T_6,size/2);
        matMul_Strassen_S<T>(T_1,T_6,Q_1,size/2);

        MatOp::AddMat<T>(A_21,A_22,T_2,size/2);
        matMul_Strassen_S<T>(T_2,B_11,Q_2,size/2);

        MatOp::SubMat<T>(B_12,B_22,T_7,size/2);
        matMul_Strassen_S<T>(A_11,T_7,Q_3,size/2);

        MatOp::SubMat<T>(B_21,B_11,T_8,size/2);
        matMul_Strassen_S<T>(A_22,T_8,Q_4,size/2);

        MatOp::AddMat<T>(A_11,A_12,T_3,size/2);
        matMul_Strassen_S<T>(T_3,B_22,Q_5,size/2);

        MatOp::SubMat<T>(A_21,A_11,T_4,size/2);
        MatOp::AddMat<T>(B_11,B_12,T_9,size/2);
        matMul_Strassen_S<T>(T_4,T_9,Q_6,size/2);

        MatOp::SubMat<T>(A_12,A_22,T_5,size/2);
        MatOp::AddMat<T>(B_11,B_22,T_10,size/2);
        matMul_Strassen_S<T>(T_5,T_10,Q_7,size/2);

        Utility::FreeMemory2D<T>(T_1);
        Utility::FreeMemory2D<T>(T_2);
        Utility::FreeMemory2D<T>(T_3);
        Utility::FreeMemory2D<T>(T_4);
        Utility::FreeMemory2D<T>(T_5);
        Utility::FreeMemory2D<T>(T_6);
        Utility::FreeMemory2D<T>(T_7);
        Utility::FreeMemory2D<T>(T_8);
        Utility::FreeMemory2D<T>(T_9);
        Utility::FreeMemory2D<T>(T_10);
        Utility::FreeMemory2D<T>(B_11);
        Utility::FreeMemory2D<T>(B_22);
        Utility::FreeMemory2D<T>(B_12);
        Utility::FreeMemory2D<T>(B_21);
        Utility::FreeMemory2D<T>(A_11);
        Utility::FreeMemory2D<T>(A_22);
        Utility::FreeMemory2D<T>(A_12);
        Utility::FreeMemory2D<T>(A_21);

        T** U_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** U_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** U_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** U_4 = Utility::AllocateMemory2D<T>(size/2,size/2);

        MatOp::AddMat<T>(Q_1,Q_4,U_1,size/2);
        MatOp::SubMat<T>(Q_5,Q_7,U_2,size/2);
        MatOp::AddMat<T>(Q_3,Q_1,U_3,size/2);
        MatOp::SubMat<T>(Q_2,Q_6,U_4,size/2);

        T** C_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** C_12 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** C_21 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** C_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
        
        MatOp::SubMat<T>(U_1,U_2,C_11,size/2);
        MatOp::AddMat<T>(Q_3,Q_5,C_12,size/2);
        MatOp::AddMat<T>(Q_2,Q_4,C_21,size/2);
        MatOp::SubMat<T>(U_3,U_4,C_22,size/2);

        Utility::FreeMemory2D<T>(Q_1);
        Utility::FreeMemory2D<T>(Q_2);
        Utility::FreeMemory2D<T>(Q_3);
        Utility::FreeMemory2D<T>(Q_4);
        Utility::FreeMemory2D<T>(Q_5);
        Utility::FreeMemory2D<T>(Q_6);
        Utility::FreeMemory2D<T>(Q_7);

        Utility::FreeMemory2D<T>(U_1);
        Utility::FreeMemory2D<T>(U_2);
        Utility::FreeMemory2D<T>(U_3);
        Utility::FreeMemory2D<T>(U_4);

        MatOp::join2D<T>(C_11,C,size/2,0,0);
        MatOp::join2D<T>(C_12,C,size/2,0,size/2);
        MatOp::join2D<T>(C_21,C,size/2,size/2,0);
        MatOp::join2D<T>(C_22,C,size/2,size/2,size/2);

        Utility::FreeMemory2D<T>(C_11);
        Utility::FreeMemory2D<T>(C_12);
        Utility::FreeMemory2D<T>(C_21);
        Utility::FreeMemory2D<T>(C_22);
    }
}

template <typename T>
static void matMul_Strassen_P_body(T**& A, T**& B, T**& C, int size)
{
    if (size <= THRESHOLD)
    {
        MatOp::matMul_Naive<T>(A, B, C, size);
    } else 
    {
        T** T_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_4 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_5 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_6 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_7 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_8 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_9 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** T_10 = Utility::AllocateMemory2D<T>(size/2,size/2);

        T** B_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** B_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** B_12 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** B_21 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** A_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** A_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** A_12 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** A_21 = Utility::AllocateMemory2D<T>(size/2,size/2);

        T** Q_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_4 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_5 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_6 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** Q_7 = Utility::AllocateMemory2D<T>(size/2,size/2);

        MatOp::split2D<T>(A,A_11,size/2,0,0);
        MatOp::split2D<T>(A,A_12,size/2,0,size/2);
        MatOp::split2D<T>(A,A_21,size/2,size/2,0);
        MatOp::split2D<T>(A,A_22,size/2,size/2,size/2);
        MatOp::split2D<T>(B,B_11,size/2,0,0);
        MatOp::split2D<T>(B,B_12,size/2,0,size/2);
        MatOp::split2D<T>(B,B_21,size/2,size/2,0);
        MatOp::split2D<T>(B,B_22,size/2,size/2,size/2);

        #pragma omp task
        {
            MatOp::MatOp::AddMat<T>(A_11,A_22,T_1,size/2);
            MatOp::AddMat<T>(B_11,B_22,T_6,size/2);
            matMul_Strassen_P_body<T>(T_1,T_6,Q_1,size/2);
        }

        #pragma omp task
        {
        MatOp::AddMat<T>(A_21,A_22,T_2,size/2);
        matMul_Strassen_P_body<T>(T_2,B_11,Q_2,size/2);
        }

        #pragma omp task
        {
        MatOp::SubMat<T>(B_12,B_22,T_7,size/2);
        matMul_Strassen_P_body<T>(A_11,T_7,Q_3,size/2);
        }

        #pragma omp task
        {
        MatOp::SubMat<T>(B_21,B_11,T_8,size/2);
        matMul_Strassen_P_body<T>(A_22,T_8,Q_4,size/2);
        }

        #pragma omp task
        {
        MatOp::AddMat<T>(A_11,A_12,T_3,size/2);
        matMul_Strassen_P_body<T>(T_3,B_22,Q_5,size/2);
        }

        #pragma omp task
        {
        MatOp::SubMat<T>(A_21,A_11,T_4,size/2);
        MatOp::AddMat<T>(B_11,B_12,T_9,size/2);
        matMul_Strassen_P_body<T>(T_4,T_9,Q_6,size/2);
        }

        #pragma omp task
        {
        MatOp::SubMat<T>(A_12,A_22,T_5,size/2);
        MatOp::AddMat<T>(B_11,B_22,T_10,size/2);
        matMul_Strassen_P_body<T>(T_5,T_10,Q_7,size/2);
        }
        #pragma omp taskwait

        Utility::FreeMemory2D<T>(T_1);
        Utility::FreeMemory2D<T>(T_2);
        Utility::FreeMemory2D<T>(T_3);
        Utility::FreeMemory2D<T>(T_4);
        Utility::FreeMemory2D<T>(T_5);
        Utility::FreeMemory2D<T>(T_6);
        Utility::FreeMemory2D<T>(T_7);
        Utility::FreeMemory2D<T>(T_8);
        Utility::FreeMemory2D<T>(T_9);
        Utility::FreeMemory2D<T>(T_10);
        Utility::FreeMemory2D<T>(B_11);
        Utility::FreeMemory2D<T>(B_22);
        Utility::FreeMemory2D<T>(B_12);
        Utility::FreeMemory2D<T>(B_21);
        Utility::FreeMemory2D<T>(A_11);
        Utility::FreeMemory2D<T>(A_22);
        Utility::FreeMemory2D<T>(A_12);
        Utility::FreeMemory2D<T>(A_21);

        T** U_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** U_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** U_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** U_4 = Utility::AllocateMemory2D<T>(size/2,size/2);

        MatOp::AddMat<T>(Q_1,Q_4,U_1,size/2);
        MatOp::SubMat<T>(Q_5,Q_7,U_2,size/2);
        MatOp::AddMat<T>(Q_3,Q_1,U_3,size/2);
        MatOp::SubMat<T>(Q_2,Q_6,U_4,size/2);

        T** C_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** C_12 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** C_21 = Utility::AllocateMemory2D<T>(size/2,size/2);
        T** C_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
        
        MatOp::SubMat<T>(U_1,U_2,C_11,size/2);
        MatOp::AddMat<T>(Q_3,Q_5,C_12,size/2);
        MatOp::AddMat<T>(Q_2,Q_4,C_21,size/2);
        MatOp::SubMat<T>(U_3,U_4,C_22,size/2);

        Utility::FreeMemory2D<T>(Q_1);
        Utility::FreeMemory2D<T>(Q_2);
        Utility::FreeMemory2D<T>(Q_3);
        Utility::FreeMemory2D<T>(Q_4);
        Utility::FreeMemory2D<T>(Q_5);
        Utility::FreeMemory2D<T>(Q_6);
        Utility::FreeMemory2D<T>(Q_7);

        Utility::FreeMemory2D<T>(U_1);
        Utility::FreeMemory2D<T>(U_2);
        Utility::FreeMemory2D<T>(U_3);
        Utility::FreeMemory2D<T>(U_4);

        MatOp::join2D<T>(C_11,C,size/2,0,0);
        MatOp::join2D<T>(C_12,C,size/2,0,size/2);
        MatOp::join2D<T>(C_21,C,size/2,size/2,0);
        MatOp::join2D<T>(C_22,C,size/2,size/2,size/2);

        Utility::FreeMemory2D<T>(C_11);
        Utility::FreeMemory2D<T>(C_12);
        Utility::FreeMemory2D<T>(C_21);
        Utility::FreeMemory2D<T>(C_22);
    }
}

template <typename T>
void matMul_Strassen_P(T**& matA, T**& matB, T**& matC, int size)
{
    #pragma omp parallel 
    {
        #pragma omp single 
        {
            matMul_Strassen_P_body(matA,matB,matC,size);
        }
    }
}

#endif