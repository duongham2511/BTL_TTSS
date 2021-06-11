#ifndef MATOP_H
#define MATOP_H

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "Utils.h"


class MatOp {
    public:

    template <typename T>
    static void matMul_Naive(T** matA, T** matB, T** matC, int rowA, int colA_rowB, int colB)
    {
        for (int i = 0; i < rowA;i++){
            for (int j = 0; j <colB; j++){
                matC[i][j] = 0;
                for (int k = 0; k <colA_rowB; k++)
                {
                    matC[i][j] += matA[i][k] * matB[k][j];
                }
            }
        }
    }

    template <typename T>
    static void pad2D (T**& in_mat, int og_row, int og_col, int new_row, int new_col)
    {
        T** temp_mat = Utility::AllocateMemory2D<T>(new_row,new_col);
        for (int i = 0; i < new_row; i++)
        {
            for (int j = 0; j < new_col; j++)
            {
                temp_mat[i][j] = ((i < og_row) && (j < og_col))? in_mat[i][j] : 0;
            }
        }
        Utility::FreeMemory2D<T>(in_mat);
        in_mat = temp_mat;
        // printMat(in_mat,new_size,new_size);
    }

    template <typename T>
    static void shrink2D(T**& in_mat, int og_row, int og_col, int new_row, int new_col)
    {
        T** temp_mat = Utility::AllocateMemory2D<T>(new_row,new_col);
        int i,j;
        for (i = 0; i <new_row; i++)
            for (j = 0; j < new_col; j++)
                temp_mat[i][j] = in_mat[i][j];
        Utility::FreeMemory2D<T>(in_mat);
        in_mat = temp_mat;
    }

    template <typename T>
    static void split2D(T** in_mat, T** out_mat,int size_out, int offset_row, int offset_col)
    {
        int i,j;
        for (i = 0; i < size_out; i++)
            for (j = 0; j < size_out; j++)
            {
                out_mat[i][j] = in_mat[i+offset_row][j+offset_col];
            }
    }

    template <typename T>
    static void join2D(T** in_mat, T** out_mat, int size_in, int offset_row, int offset_col)
    {
        int i,j;
        for (i = 0; i <size_in; i++)
            for (j = 0; j <size_in; j++)
            {
                out_mat[i+offset_row][j+offset_col] = in_mat[i][j];
            }
    }

    template <typename T>
    static void AddMat(T** matA, T** matB, T**mat_out, int size)
    {
        int i,j;
        for (i = 0; i <size; i++)
            for (j = 0; j <size; j++)
            {
                mat_out[i][j] = matA[i][j] + matB[i][j];
            }
    }

    template <typename T>
    static void SubMat(T** matA, T** matB, T**mat_out, int size)
    {
        int i,j;
        for (i = 0; i <size; i++)
            for (j = 0; j <size; j++)
            {
                mat_out[i][j] = matA[i][j] - matB[i][j];
            }
    }
};

#endif