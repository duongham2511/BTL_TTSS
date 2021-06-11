#ifndef STRASSEN_H
#define STRASSEN_H

#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "MatOp.h"
#include "Utils.h"

// template <typename T>
// void matMul_Strassen_parallel(T** matA, T** matB, T** matC, int size)
// {
//     if (size <= 64)
//     {
//         MatOp::matMul_Naive<T>(matA,matB,matC,size);
//     } else{
//         T** T_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_4 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_5 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_6 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_7 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_8 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_9 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_10 = Utility::AllocateMemory2D<T>(size/2,size/2);

//         T** B_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** B_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** A_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** A_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         for (int i = 0; i < size/2; i++){
//             for (int j = 0; j < size/2; j++){
//                 B_11[i][j] = 0;
//                 B_22[i][j] = 0;
//                 A_11[i][j] = 0;
//                 A_22[i][j] = 0;
//             }
//         }

//         T** Q_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_4 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_5 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_6 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_7 = Utility::AllocateMemory2D<T>(size/2,size/2);

//         #pragma omp task
//         {
//             MatOp::AddMatSection(matA,0,0,matA,size/2,size/2,T_1,0,0,size/2,size/2);
//             //printf("\nT_1 = \n");
//             //printMat(T_1,size/2,size/2);
//             MatOp::AddMatSection(matB,0,0,matB,size/2,size/2,T_6,0,0,size/2,size/2);
//             //printf("\nT_6 = \n");
//             //printMat(T_6,size/2,size/2);
//             matMul_Strassen_parallel(T_1,T_6,Q_1,size/2);
//         }
//         #pragma omp task
//         {
//             MatOp::AddMatSection(matA,size/2,0,matA,size/2,size/2,T_2,0,0,size/2,size/2);
//             //printf("\nT_2 = \n");
//             //printMat(T_2,size/2,size/2);
//             MatOp::AddMatSection(matB,0,0,B_11,0,0,B_11,0,0,size/2,size/2);
//             //printMat(B_11,size/2,size/2);
//             matMul_Strassen_parallel(T_2,B_11,Q_2,size/2);
//         }
//         #pragma omp task
//         {
//             MatOp::AddMatSection(matA,0,0,A_11,0,0,A_11,0,0,size/2,size/2);
//             //printMat(A_11,size/2,size/2);
//             MatOp::SubMatSection(matB,0,size/2,matB,size/2,size/2,T_7,0,0,size/2,size/2);
//             //printf("\nT_7 = \n");
//             //printMat(T_7,size/2,size/2);
//             matMul_Strassen_parallel(A_11,T_7,Q_3,size/2);
//         }
//         #pragma omp task
//         {
//             MatOp::AddMatSection(matA,size/2,size/2,A_22,0,0,A_22,0,0,size/2,size/2);
//             //printMat(A_22,size/2,size/2);
//             MatOp::SubMatSection(matB,size/2,0,matB,0,0,T_8,0,0,size/2,size/2);
//             //printf("\nT_8 = \n");
//             //printMat(T_8,size/2,size/2);
//             matMul_Strassen_parallel(A_22,T_8,Q_4,size/2);
//         }
//         #pragma omp task
//         {
//             MatOp::AddMatSection(matA,0,0,matA,0,size/2,T_3,0,0,size/2,size/2);
//             //printf("\nT_3 = \n");
//             //printMat(T_3,size/2,size/2);
//             MatOp::AddMatSection(matB,size/2,size/2,B_22,0,0,B_22,0,0,size/2,size/2);
//             //printMat(B_22,size/2,size/2);
//             matMul_Strassen_parallel(T_3,B_22,Q_5,size/2);
//         }
//         #pragma omp task
//         {
//             MatOp::SubMatSection(matA,size/2,0,matA,0,0,T_4,0,0,size/2,size/2);
//             //printf("\nT_4 = \n");
//             //printMat(T_4,size/2,size/2);
//             MatOp::AddMatSection(matB,0,0,matB,0,size/2,T_9,0,0,size/2,size/2);
//             //printf("\nT_9 = \n");
//             //printMat(T_9,size/2,size/2);
//             matMul_Strassen_parallel(T_4,T_9,Q_6,size/2);
//         }
//         #pragma omp task
//         {
//             MatOp::SubMatSection(matA,0,size/2,matA,size/2,size/2,T_5,0,0,size/2,size/2);
//             //printf("\nT_5 = \n");
//             //printMat(T_5,size/2,size/2);
//             MatOp::AddMatSection(matB,size/2,0,matB,size/2,size/2,T_10,0,0,size/2,size/2);
//             //printf("\nT_10 = \n");
//             //printMat(T_10,size/2,size/2);
//             matMul_Strassen_parallel(T_5,T_10,Q_7,size/2);
//         }
//         #pragma omp taskwait

//         Utility::FreeMemory2D<T>(T_1);
//         Utility::FreeMemory2D<T>(T_2);
//         Utility::FreeMemory2D<T>(T_3);
//         Utility::FreeMemory2D<T>(T_4);
//         Utility::FreeMemory2D<T>(T_5);
//         Utility::FreeMemory2D<T>(T_6);
//         Utility::FreeMemory2D<T>(T_7);
//         Utility::FreeMemory2D<T>(T_8);
//         Utility::FreeMemory2D<T>(T_9);
//         Utility::FreeMemory2D<T>(T_10);
//         Utility::FreeMemory2D<T>(B_11);
//         Utility::FreeMemory2D<T>(B_22);
//         Utility::FreeMemory2D<T>(A_11);
//         Utility::FreeMemory2D<T>(A_22);

//         T** U_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** U_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** U_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** U_4 = Utility::AllocateMemory2D<T>(size/2,size/2);

//         MatOp::AddMatSection(Q_1,0,0,Q_4,0,0,U_1,0,0,size/2,size/2);
//         MatOp::SubMatSection(Q_5,0,0,Q_7,0,0,U_2,0,0,size/2,size/2);
//         MatOp::AddMatSection(Q_3,0,0,Q_1,0,0,U_3,0,0,size/2,size/2);
//         MatOp::SubMatSection(Q_2,0,0,Q_6,0,0,U_4,0,0,size/2,size/2);

//         MatOp::SubMatSection(U_1,0,0,U_2,0,0,matC,0,0,size/2,size/2);
//         MatOp::AddMatSection(Q_3,0,0,Q_5,0,0,matC,0,size/2,size/2,size/2);
//         MatOp::AddMatSection(Q_2,0,0,Q_4,0,0,matC,size/2,0,size/2,size/2);
//         MatOp::SubMatSection(U_3,0,0,U_4,0,0,matC,size/2,size/2,size/2,size/2);

//         Utility::FreeMemory2D<T>(Q_1);
//         Utility::FreeMemory2D<T>(Q_2);
//         Utility::FreeMemory2D<T>(Q_3);
//         Utility::FreeMemory2D<T>(Q_4);
//         Utility::FreeMemory2D<T>(Q_5);
//         Utility::FreeMemory2D<T>(Q_6);
//         Utility::FreeMemory2D<T>(Q_7);

//         Utility::FreeMemory2D<T>(U_1);
//         Utility::FreeMemory2D<T>(U_2);
//         Utility::FreeMemory2D<T>(U_3);
//         Utility::FreeMemory2D<T>(U_4);
//     }
// }

// template <typename T>
// void matMul_Strassen_v2(T** matA, T** matB, T** matC, int size)
// {
//     if (size <= 32)
//     {
//         MatOp::matMul_Naive<T>(matA,matB,matC,size);
//     } else{
//         T** T_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_4 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_5 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_6 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_7 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_8 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_9 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_10 = Utility::AllocateMemory2D<T>(size/2,size/2);

//         T** B_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** B_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** A_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** A_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         for (int i = 0; i < size/2; i++){
//             for (int j = 0; j < size/2; j++){
//                 B_11[i][j] = 0;
//                 B_22[i][j] = 0;
//                 A_11[i][j] = 0;
//                 A_22[i][j] = 0;
//             }
//         }

//         MatOp::AddMatSection(matA,0,0,matA,size/2,size/2,T_1,0,0,size/2,size/2);
//         //printf("\nT_1 = \n");
//         //printMat(T_1,size/2,size/2);
//         MatOp::AddMatSection(matA,size/2,0,matA,size/2,size/2,T_2,0,0,size/2,size/2);
//         //printf("\nT_2 = \n");
//         //printMat(T_2,size/2,size/2);
//         MatOp::AddMatSection(matA,0,0,matA,0,size/2,T_3,0,0,size/2,size/2);
//         //printf("\nT_3 = \n");
//         //printMat(T_3,size/2,size/2);
//         MatOp::SubMatSection(matA,size/2,0,matA,0,0,T_4,0,0,size/2,size/2);
//         //printf("\nT_4 = \n");
//         //printMat(T_4,size/2,size/2);
//         MatOp::SubMatSection(matA,0,size/2,matA,size/2,size/2,T_5,0,0,size/2,size/2);
//         //printf("\nT_5 = \n");
//         //printMat(T_5,size/2,size/2);
//         MatOp::AddMatSection(matB,0,0,matB,size/2,size/2,T_6,0,0,size/2,size/2);
//         //printf("\nT_6 = \n");
//         //printMat(T_6,size/2,size/2);
//         MatOp::SubMatSection(matB,0,size/2,matB,size/2,size/2,T_7,0,0,size/2,size/2);
//         //printf("\nT_7 = \n");
//         //printMat(T_7,size/2,size/2);
//         MatOp::SubMatSection(matB,size/2,0,matB,0,0,T_8,0,0,size/2,size/2);
//         //printf("\nT_8 = \n");
//         //printMat(T_8,size/2,size/2);
//         MatOp::AddMatSection(matB,0,0,matB,0,size/2,T_9,0,0,size/2,size/2);
//         //printf("\nT_9 = \n");
//         //printMat(T_9,size/2,size/2);
//         MatOp::AddMatSection(matB,size/2,0,matB,size/2,size/2,T_10,0,0,size/2,size/2);
//         //printf("\nT_10 = \n");
//         //printMat(T_10,size/2,size/2);

//         MatOp::AddMatSection(matB,0,0,B_11,0,0,B_11,0,0,size/2,size/2);
//         //printMat(B_11,size/2,size/2);
//         MatOp::AddMatSection(matB,size/2,size/2,B_22,0,0,B_22,0,0,size/2,size/2);
//         //printMat(B_22,size/2,size/2);
//         MatOp::AddMatSection(matA,0,0,A_11,0,0,A_11,0,0,size/2,size/2);
//         //printMat(A_11,size/2,size/2);
//         MatOp::AddMatSection(matA,size/2,size/2,A_22,0,0,A_22,0,0,size/2,size/2);
//         //printMat(A_22,size/2,size/2);

//         T** Q_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_4 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_5 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_6 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_7 = Utility::AllocateMemory2D<T>(size/2,size/2);

//         matMul_Strassen_v2(T_1,T_6,Q_1,size/2);
//         matMul_Strassen_v2(T_2,B_11,Q_2,size/2);
//         matMul_Strassen_v2(A_11,T_7,Q_3,size/2);
//         matMul_Strassen_v2(A_22,T_8,Q_4,size/2);
//         matMul_Strassen_v2(T_3,B_22,Q_5,size/2);
//         matMul_Strassen_v2(T_4,T_9,Q_6,size/2);
//         matMul_Strassen_v2(T_5,T_10,Q_7,size/2);

//         Utility::FreeMemory2D<T>(T_1);
//         Utility::FreeMemory2D<T>(T_2);
//         Utility::FreeMemory2D<T>(T_3);
//         Utility::FreeMemory2D<T>(T_4);
//         Utility::FreeMemory2D<T>(T_5);
//         Utility::FreeMemory2D<T>(T_6);
//         Utility::FreeMemory2D<T>(T_7);
//         Utility::FreeMemory2D<T>(T_8);
//         Utility::FreeMemory2D<T>(T_9);
//         Utility::FreeMemory2D<T>(T_10);
//         Utility::FreeMemory2D<T>(B_11);
//         Utility::FreeMemory2D<T>(B_22);
//         Utility::FreeMemory2D<T>(A_11);
//         Utility::FreeMemory2D<T>(A_22);

//         T** U_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** U_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** U_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** U_4 = Utility::AllocateMemory2D<T>(size/2,size/2);

//         MatOp::AddMatSection(Q_1,0,0,Q_4,0,0,U_1,0,0,size/2,size/2);
//         MatOp::SubMatSection(Q_5,0,0,Q_7,0,0,U_2,0,0,size/2,size/2);
//         MatOp::AddMatSection(Q_3,0,0,Q_1,0,0,U_3,0,0,size/2,size/2);
//         MatOp::SubMatSection(Q_2,0,0,Q_6,0,0,U_4,0,0,size/2,size/2);

//         MatOp::SubMatSection(U_1,0,0,U_2,0,0,matC,0,0,size/2,size/2);
//         MatOp::AddMatSection(Q_3,0,0,Q_5,0,0,matC,0,size/2,size/2,size/2);
//         MatOp::AddMatSection(Q_2,0,0,Q_4,0,0,matC,size/2,0,size/2,size/2);
//         MatOp::SubMatSection(U_3,0,0,U_4,0,0,matC,size/2,size/2,size/2,size/2);

//         Utility::FreeMemory2D<T>(Q_1);
//         Utility::FreeMemory2D<T>(Q_2);
//         Utility::FreeMemory2D<T>(Q_3);
//         Utility::FreeMemory2D<T>(Q_4);
//         Utility::FreeMemory2D<T>(Q_5);
//         Utility::FreeMemory2D<T>(Q_6);
//         Utility::FreeMemory2D<T>(Q_7);

//         Utility::FreeMemory2D<T>(U_1);
//         Utility::FreeMemory2D<T>(U_2);
//         Utility::FreeMemory2D<T>(U_3);
//         Utility::FreeMemory2D<T>(U_4);
//     }
// }

// template <typename T>
// void matMul_Strassen_serial(T** matA, T** matB, T** matC, int size)
// {
//     if (size <= 64)
//     {
//         MatOp::matMul_Naive<T>(matA,matB,matC,size);
//     } else{
//         T** T_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_4 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_5 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_6 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_7 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_8 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_9 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** T_10 = Utility::AllocateMemory2D<T>(size/2,size/2);

//         T** B_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** B_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** A_11 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** A_22 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         for (int i = 0; i < size/2; i++){
//             for (int j = 0; j < size/2; j++){
//                 B_11[i][j] = 0;
//                 B_22[i][j] = 0;
//                 A_11[i][j] = 0;
//                 A_22[i][j] = 0;
//             }
//         }

//         MatOp::AddMatSection(matA,0,0,matA,size/2,size/2,T_1,0,0,size/2,size/2);
//         //printf("\nT_1 = \n");
//         //printMat(T_1,size/2,size/2);
//         MatOp::AddMatSection(matA,size/2,0,matA,size/2,size/2,T_2,0,0,size/2,size/2);
//         //printf("\nT_2 = \n");
//         //printMat(T_2,size/2,size/2);
//         MatOp::AddMatSection(matA,0,0,matA,0,size/2,T_3,0,0,size/2,size/2);
//         //printf("\nT_3 = \n");
//         //printMat(T_3,size/2,size/2);
//         MatOp::SubMatSection(matA,size/2,0,matA,0,0,T_4,0,0,size/2,size/2);
//         //printf("\nT_4 = \n");
//         //printMat(T_4,size/2,size/2);
//         MatOp::SubMatSection(matA,0,size/2,matA,size/2,size/2,T_5,0,0,size/2,size/2);
//         //printf("\nT_5 = \n");
//         //printMat(T_5,size/2,size/2);
//         MatOp::AddMatSection(matB,0,0,matB,size/2,size/2,T_6,0,0,size/2,size/2);
//         //printf("\nT_6 = \n");
//         //printMat(T_6,size/2,size/2);
//         MatOp::SubMatSection(matB,0,size/2,matB,size/2,size/2,T_7,0,0,size/2,size/2);
//         //printf("\nT_7 = \n");
//         //printMat(T_7,size/2,size/2);
//         MatOp::SubMatSection(matB,size/2,0,matB,0,0,T_8,0,0,size/2,size/2);
//         //printf("\nT_8 = \n");
//         //printMat(T_8,size/2,size/2);
//         MatOp::AddMatSection(matB,0,0,matB,0,size/2,T_9,0,0,size/2,size/2);
//         //printf("\nT_9 = \n");
//         //printMat(T_9,size/2,size/2);
//         MatOp::AddMatSection(matB,size/2,0,matB,size/2,size/2,T_10,0,0,size/2,size/2);
//         //printf("\nT_10 = \n");
//         //printMat(T_10,size/2,size/2);

//         MatOp::AddMatSection(matB,0,0,B_11,0,0,B_11,0,0,size/2,size/2);
//         //printMat(B_11,size/2,size/2);
//         MatOp::AddMatSection(matB,size/2,size/2,B_22,0,0,B_22,0,0,size/2,size/2);
//         //printMat(B_22,size/2,size/2);
//         MatOp::AddMatSection(matA,0,0,A_11,0,0,A_11,0,0,size/2,size/2);
//         //printMat(A_11,size/2,size/2);
//         MatOp::AddMatSection(matA,size/2,size/2,A_22,0,0,A_22,0,0,size/2,size/2);
//         //printMat(A_22,size/2,size/2);

//         T** Q_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_4 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_5 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_6 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** Q_7 = Utility::AllocateMemory2D<T>(size/2,size/2);

//         matMul_Strassen_serial(T_1,T_6,Q_1,size/2);
//         matMul_Strassen_serial(T_2,B_11,Q_2,size/2);
//         matMul_Strassen_serial(A_11,T_7,Q_3,size/2);
//         matMul_Strassen_serial(A_22,T_8,Q_4,size/2);
//         matMul_Strassen_serial(T_3,B_22,Q_5,size/2);
//         matMul_Strassen_serial(T_4,T_9,Q_6,size/2);
//         matMul_Strassen_serial(T_5,T_10,Q_7,size/2);

//         Utility::FreeMemory2D<T>(T_1);
//         Utility::FreeMemory2D<T>(T_2);
//         Utility::FreeMemory2D<T>(T_3);
//         Utility::FreeMemory2D<T>(T_4);
//         Utility::FreeMemory2D<T>(T_5);
//         Utility::FreeMemory2D<T>(T_6);
//         Utility::FreeMemory2D<T>(T_7);
//         Utility::FreeMemory2D<T>(T_8);
//         Utility::FreeMemory2D<T>(T_9);
//         Utility::FreeMemory2D<T>(T_10);
//         Utility::FreeMemory2D<T>(B_11);
//         Utility::FreeMemory2D<T>(B_22);
//         Utility::FreeMemory2D<T>(A_11);
//         Utility::FreeMemory2D<T>(A_22);

//         T** U_1 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** U_2 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** U_3 = Utility::AllocateMemory2D<T>(size/2,size/2);
//         T** U_4 = Utility::AllocateMemory2D<T>(size/2,size/2);

//         MatOp::AddMatSection(Q_1,0,0,Q_4,0,0,U_1,0,0,size/2,size/2);
//         MatOp::SubMatSection(Q_5,0,0,Q_7,0,0,U_2,0,0,size/2,size/2);
//         MatOp::AddMatSection(Q_3,0,0,Q_1,0,0,U_3,0,0,size/2,size/2);
//         MatOp::SubMatSection(Q_2,0,0,Q_6,0,0,U_4,0,0,size/2,size/2);

//         MatOp::SubMatSection(U_1,0,0,U_2,0,0,matC,0,0,size/2,size/2);
//         MatOp::AddMatSection(Q_3,0,0,Q_5,0,0,matC,0,size/2,size/2,size/2);
//         MatOp::AddMatSection(Q_2,0,0,Q_4,0,0,matC,size/2,0,size/2,size/2);
//         MatOp::SubMatSection(U_3,0,0,U_4,0,0,matC,size/2,size/2,size/2,size/2);

//         Utility::FreeMemory2D<T>(Q_1);
//         Utility::FreeMemory2D<T>(Q_2);
//         Utility::FreeMemory2D<T>(Q_3);
//         Utility::FreeMemory2D<T>(Q_4);
//         Utility::FreeMemory2D<T>(Q_5);
//         Utility::FreeMemory2D<T>(Q_6);
//         Utility::FreeMemory2D<T>(Q_7);

//         Utility::FreeMemory2D<T>(U_1);
//         Utility::FreeMemory2D<T>(U_2);
//         Utility::FreeMemory2D<T>(U_3);
//         Utility::FreeMemory2D<T>(U_4);
//     }
// }

// template <typename T>
// void matMul_parallel(T**& matA, T**& matB, T**& matC, int size)
// {
//     #pragma omp parallel 
//     {
//         #pragma omp single 
//         {
//             omp_set_num_threads(8);
//             matMul_Strassen_parallel(matA,matB,matC,size);
//         }
//     }
// }

template <typename T>
void matMul_Strassen_body(T**& A, T**& B, T**& C, int size)
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
        matMul_Strassen_body<T>(T_1,T_6,Q_1,size/2);

        MatOp::AddMat<T>(A_21,A_22,T_2,size/2);
        matMul_Strassen_body<T>(T_2,B_11,Q_2,size/2);

        MatOp::SubMat<T>(B_12,B_22,T_7,size/2);
        matMul_Strassen_body<T>(A_11,T_7,Q_3,size/2);

        MatOp::SubMat<T>(B_21,B_11,T_8,size/2);
        matMul_Strassen_body<T>(A_22,T_8,Q_4,size/2);

        MatOp::AddMat<T>(A_11,A_12,T_3,size/2);
        matMul_Strassen_body<T>(T_3,B_22,Q_5,size/2);

        MatOp::SubMat<T>(A_21,A_11,T_4,size/2);
        MatOp::AddMat<T>(B_11,B_12,T_9,size/2);
        matMul_Strassen_body<T>(T_4,T_9,Q_6,size/2);

        MatOp::SubMat<T>(A_12,A_22,T_5,size/2);
        MatOp::AddMat<T>(B_11,B_22,T_10,size/2);
        matMul_Strassen_body<T>(T_5,T_10,Q_7,size/2);

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

#endif