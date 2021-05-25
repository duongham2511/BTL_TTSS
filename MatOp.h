#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

template <typename T>
void AddMatSection(T** matA, int offset_rowA, int offset_colA, T** matB, int offset_rowB, int offset_colB, T** matC, int offset_rowC, int offset_colC, int size_row, int size_col);

template <typename T>
void SubMatSection(T** matA, int offset_rowA, int offset_colA, T** matB, int offset_rowB, int offset_colB, T** matC, int offset_rowC, int offset_colC, int size_row, int size_col);

template <typename T>
void matMul_Strassen(T** matA, T** matB, T** matC, int size);

template <typename T>
void matMul_Strassen_v2(T** matA, T** matB, T** matC, int size);

template <typename T>
void matMul_Naive(T** matA, T** matB, T** matC, int size);

template <typename T>
T** AllocateMemory2D (int height, int width);

template <typename T>
void FreeMemory2D(T** arr);

int getPower2(int value);

template <typename T>
int pad2D (T**& in_mat, int og_size);

template <typename T>
void shrink2D(T**& in_mat, int og_size, int new_size);

template <typename T>
void split2D(T** in_mat, T** out_mat,int size_out, int offset_row, int offset_col);

template <typename T>
void join2D(T** in_mat, T** out_mat, int size_in, int offset_row, int offset_col);

template <typename T>
void AddMat(T** matA, T** matB, T**mat_out, int size);

template <typename T>
void SubMat(T** matA, T** matB, T**mat_out, int size);

int getPower2(int value)
{
    int pow = 1;
    while (pow < value) pow = pow*2;
    return pow;
}

void printMat(int** mat, int height, int width)
{
    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++)
        {
            printf("%3d",mat[i][j]);
            if (j == width -1) printf("\n");
            else printf(" ");
        }
    }
}

void printMat(double** mat, int height, int width)
{
    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++)
        {
            printf("%6.2f",mat[i][j]);
            if (j == width -1) printf("\n");
            else printf(" ");
        }
    }
}

template <typename T>
int pad2D (T**& in_mat, int og_size)
{
    int new_size = getPower2(og_size);
    T** temp_mat = AllocateMemory2D<T>(new_size,new_size);
    for (int i = 0; i < new_size; i++)
    {
        for (int j = 0; j < new_size; j++)
        {
            temp_mat[i][j] = ((i < og_size) && (j < og_size))? in_mat[i][j] : 0;
        }
    }
    FreeMemory2D<T>(in_mat);
    in_mat = temp_mat;
    // printMat(in_mat,new_size,new_size);
    return new_size;
}

template <typename T>
void shrink2D(T**& in_mat, int og_size, int new_size)
{
    T** temp_mat = AllocateMemory2D<T>(new_size,new_size);
    int i,j;
    for (i = 0; i <new_size; i++)
        for (j = 0; j < new_size; j++)
            temp_mat[i][j] = in_mat[i][j];
    FreeMemory2D<T>(in_mat);
    in_mat = temp_mat;
}

template <typename T>
void split2D(T** in_mat, T** out_mat,int size_out, int offset_row, int offset_col)
{
    int i,j;
    for (i = 0; i < size_out; i++)
        for (j = 0; j < size_out; j++)
        {
            out_mat[i][j] = in_mat[i+offset_row][j+offset_col];
        }
}

template <typename T>
void join2D(T** in_mat, T** out_mat, int size_in, int offset_row, int offset_col)
{
    int i,j;
    for (i = 0; i <size_in; i++)
        for (j = 0; j <size_in; j++)
        {
            out_mat[i+offset_row][j+offset_col] = in_mat[i][j];
        }
}

template <typename T>
T** AllocateMemory2D (int height, int width)
{
    T** ret = new T*[height];
    T* arr = new T[height*width];
    for (int i = 0; i < height; i++)
    {
        *(ret+i) = arr;
        arr += width;
    }
    return ret;
}

template <typename T>
void FreeMemory2D(T** arr)
{
    delete[] *arr;
    delete[] arr;
}

template <typename T>
void AddMat(T** matA, T** matB, T**mat_out, int size)
{
    int i,j;
    for (i = 0; i <size; i++)
        for (j = 0; j <size; j++)
        {
            mat_out[i][j] = matA[i][j] + matB[i][j];
        }
}

template <typename T>
void SubMat(T** matA, T** matB, T**mat_out, int size)
{
    int i,j;
    for (i = 0; i <size; i++)
        for (j = 0; j <size; j++)
        {
            mat_out[i][j] = matA[i][j] - matB[i][j];
        }
}

template <typename T>
void AddMatSection(T** matA, int offset_rowA, int offset_colA, T** matB, int offset_rowB, int offset_colB, T** matC, int offset_rowC, int offset_colC, int size_row, int size_col)
/*
    perform addition with a submatrix of matA and a submatrix of matB and store the result in submatrix of matC
*/
{
    try
    {
        if ((offset_rowA <0) || (offset_colA <0) || (offset_rowB <0) || (offset_colB <0) || (offset_rowC <0) || (offset_colC <0)) throw "Offset_lower than 0";
    }
    catch(const char* e)
    {
        std::cerr << e << '\n';
        exit(-1);
    }
    T** matBuf = AllocateMemory2D<T>(size_row,size_col);
    int i,j;
    #pragma omp parallel for default(shared) private (i,j)
    for (i = 0; i < size_row; i++)
    {
        for (j = 0; j < size_col; j++)
        {
            matBuf[i][j] = matA[i+offset_rowA][j+offset_colA] + matB[i+offset_rowB][j+offset_colB];
        }
    }
    #pragma omp parallel for default(shared) private (i,j)
    for (i = 0; i < size_row; i++)
    {
        for (j = 0; j < size_col; j++)
        {
            matC[i + offset_rowC][j + offset_colC] = matBuf[i][j];
        }
    }
    FreeMemory2D<T>(matBuf);
}

template <typename T>
void SubMatSection(T** matA, int offset_rowA, int offset_colA, T** matB, int offset_rowB, int offset_colB, T** matC, int offset_rowC, int offset_colC, int size_row, int size_col)
/*
    perform subtraction with a submatrix of matA and a submatrix of matB and store the result in submatrix of matC
*/
{
    try
    {
        if ((offset_rowA <0) || (offset_colA <0) || (offset_rowB <0) || (offset_colB <0) || (offset_rowC <0) || (offset_colC <0)) throw "Offset_lower than 0";
    }
    catch(const char* e)
    {
        std::cerr << e << '\n';
        exit(-1);
    }
    T** matBuf = AllocateMemory2D<T>(size_row,size_col);
    int i,j;
    #pragma omp parallel for default(shared) private (i,j)
    for (i = 0; i < size_row; i++)
    {
        for (j = 0; j < size_col; j++)
        {
            matBuf[i][j] = matA[i+offset_rowA][j+offset_colA] - matB[i+offset_rowB][j+offset_colB];
        }
    }
    #pragma omp parallel for default(shared) private (i,j)
    for (i = 0; i < size_row; i++)
    {
        for (j = 0; j < size_col; j++)
        {
            matC[i + offset_rowC][j + offset_colC] = matBuf[i][j];
        }
    }
    FreeMemory2D<T>(matBuf);
}

template <typename T>
void matMul_Strassen(T** matA, T** matB, T** matC, int size)
{
    if (size ==1)
    {
        matC[0][0] = matA[0][0] * matB[0][0];
    } else{
        T** T_1 = AllocateMemory2D<T>(size/2,size/2);
        T** T_2 = AllocateMemory2D<T>(size/2,size/2);
        T** T_3 = AllocateMemory2D<T>(size/2,size/2);
        T** T_4 = AllocateMemory2D<T>(size/2,size/2);
        T** T_5 = AllocateMemory2D<T>(size/2,size/2);
        T** T_6 = AllocateMemory2D<T>(size/2,size/2);
        T** T_7 = AllocateMemory2D<T>(size/2,size/2);
        T** T_8 = AllocateMemory2D<T>(size/2,size/2);
        T** T_9 = AllocateMemory2D<T>(size/2,size/2);
        T** T_10 = AllocateMemory2D<T>(size/2,size/2);

        T** B_11 = AllocateMemory2D<T>(size/2,size/2);
        T** B_22 = AllocateMemory2D<T>(size/2,size/2);
        T** A_11 = AllocateMemory2D<T>(size/2,size/2);
        T** A_22 = AllocateMemory2D<T>(size/2,size/2);
        for (int i = 0; i < size/2; i++){
            for (int j = 0; j < size/2; j++){
                B_11[i][j] = 0;
                B_22[i][j] = 0;
                A_11[i][j] = 0;
                A_22[i][j] = 0;
            }
        }

        AddMatSection(matA,0,0,matA,size/2,size/2,T_1,0,0,size/2,size/2);
        //printf("\nT_1 = \n");
        //printMat(T_1,size/2,size/2);
        AddMatSection(matA,size/2,0,matA,size/2,size/2,T_2,0,0,size/2,size/2);
        //printf("\nT_2 = \n");
        //printMat(T_2,size/2,size/2);
        AddMatSection(matA,0,0,matA,0,size/2,T_3,0,0,size/2,size/2);
        //printf("\nT_3 = \n");
        //printMat(T_3,size/2,size/2);
        SubMatSection(matA,size/2,0,matA,0,0,T_4,0,0,size/2,size/2);
        //printf("\nT_4 = \n");
        //printMat(T_4,size/2,size/2);
        SubMatSection(matA,0,size/2,matA,size/2,size/2,T_5,0,0,size/2,size/2);
        //printf("\nT_5 = \n");
        //printMat(T_5,size/2,size/2);
        AddMatSection(matB,0,0,matB,size/2,size/2,T_6,0,0,size/2,size/2);
        //printf("\nT_6 = \n");
        //printMat(T_6,size/2,size/2);
        SubMatSection(matB,0,size/2,matB,size/2,size/2,T_7,0,0,size/2,size/2);
        //printf("\nT_7 = \n");
        //printMat(T_7,size/2,size/2);
        SubMatSection(matB,size/2,0,matB,0,0,T_8,0,0,size/2,size/2);
        //printf("\nT_8 = \n");
        //printMat(T_8,size/2,size/2);
        AddMatSection(matB,0,0,matB,0,size/2,T_9,0,0,size/2,size/2);
        //printf("\nT_9 = \n");
        //printMat(T_9,size/2,size/2);
        AddMatSection(matB,size/2,0,matB,size/2,size/2,T_10,0,0,size/2,size/2);
        //printf("\nT_10 = \n");
        //printMat(T_10,size/2,size/2);

        AddMatSection(matB,0,0,B_11,0,0,B_11,0,0,size/2,size/2);
        //printMat(B_11,size/2,size/2);
        AddMatSection(matB,size/2,size/2,B_22,0,0,B_22,0,0,size/2,size/2);
        //printMat(B_22,size/2,size/2);
        AddMatSection(matA,0,0,A_11,0,0,A_11,0,0,size/2,size/2);
        //printMat(A_11,size/2,size/2);
        AddMatSection(matA,size/2,size/2,A_22,0,0,A_22,0,0,size/2,size/2);
        //printMat(A_22,size/2,size/2);

        T** Q_1 = AllocateMemory2D<T>(size/2,size/2);
        T** Q_2 = AllocateMemory2D<T>(size/2,size/2);
        T** Q_3 = AllocateMemory2D<T>(size/2,size/2);
        T** Q_4 = AllocateMemory2D<T>(size/2,size/2);
        T** Q_5 = AllocateMemory2D<T>(size/2,size/2);
        T** Q_6 = AllocateMemory2D<T>(size/2,size/2);
        T** Q_7 = AllocateMemory2D<T>(size/2,size/2);

        matMul_Strassen(T_1,T_6,Q_1,size/2);
        matMul_Strassen(T_2,B_11,Q_2,size/2);
        matMul_Strassen(A_11,T_7,Q_3,size/2);
        matMul_Strassen(A_22,T_8,Q_4,size/2);
        matMul_Strassen(T_3,B_22,Q_5,size/2);
        matMul_Strassen(T_4,T_9,Q_6,size/2);
        matMul_Strassen(T_5,T_10,Q_7,size/2);

        FreeMemory2D<T>(T_1);
        FreeMemory2D<T>(T_2);
        FreeMemory2D<T>(T_3);
        FreeMemory2D<T>(T_4);
        FreeMemory2D<T>(T_5);
        FreeMemory2D<T>(T_6);
        FreeMemory2D<T>(T_7);
        FreeMemory2D<T>(T_8);
        FreeMemory2D<T>(T_9);
        FreeMemory2D<T>(T_10);
        FreeMemory2D<T>(B_11);
        FreeMemory2D<T>(B_22);
        FreeMemory2D<T>(A_11);
        FreeMemory2D<T>(A_22);

        T** U_1 = AllocateMemory2D<T>(size/2,size/2);
        T** U_2 = AllocateMemory2D<T>(size/2,size/2);
        T** U_3 = AllocateMemory2D<T>(size/2,size/2);
        T** U_4 = AllocateMemory2D<T>(size/2,size/2);

        AddMatSection(Q_1,0,0,Q_4,0,0,U_1,0,0,size/2,size/2);
        SubMatSection(Q_5,0,0,Q_7,0,0,U_2,0,0,size/2,size/2);
        AddMatSection(Q_3,0,0,Q_1,0,0,U_3,0,0,size/2,size/2);
        SubMatSection(Q_2,0,0,Q_6,0,0,U_4,0,0,size/2,size/2);

        SubMatSection(U_1,0,0,U_2,0,0,matC,0,0,size/2,size/2);
        AddMatSection(Q_3,0,0,Q_5,0,0,matC,0,size/2,size/2,size/2);
        AddMatSection(Q_2,0,0,Q_4,0,0,matC,size/2,0,size/2,size/2);
        SubMatSection(U_3,0,0,U_4,0,0,matC,size/2,size/2,size/2,size/2);

        FreeMemory2D<T>(Q_1);
        FreeMemory2D<T>(Q_2);
        FreeMemory2D<T>(Q_3);
        FreeMemory2D<T>(Q_4);
        FreeMemory2D<T>(Q_5);
        FreeMemory2D<T>(Q_6);
        FreeMemory2D<T>(Q_7);

        FreeMemory2D<T>(U_1);
        FreeMemory2D<T>(U_2);
        FreeMemory2D<T>(U_3);
        FreeMemory2D<T>(U_4);
    }
}

template <typename T>
void matMul_Naive(T** matA, T** matB, T** matC, int size)
{
    for (int i = 0; i <size;i++){
        for (int j = 0; j <size; j++){
            matC[i][j] = 0;
            for (int k = 0; k <size; k++)
            {
                matC[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }
}

template <typename T>
void matMul_Strassen_v2(T** matA, T** matB, T** matC, int size)
{
    // NOTE: INCOMPLETE
    if (size ==1)
    {
        matC[0][0] = matA[0][0] * matB[0][0];
    } else{
        T** X11 = AllocateMemory2D<T>(size/2, size/2);
        T** X12 = AllocateMemory2D<T>(size/2, size/2);
        T** X21 = AllocateMemory2D<T>(size/2, size/2);
        T** X22 = AllocateMemory2D<T>(size/2, size/2);

        T** Y11 = AllocateMemory2D<T>(size/2, size/2);
        T** Y12 = AllocateMemory2D<T>(size/2, size/2);
        T** Y21 = AllocateMemory2D<T>(size/2, size/2);
        T** Y22 = AllocateMemory2D<T>(size/2, size/2);

        T** Z11 = AllocateMemory2D<T>(size/2, size/2);
        T** Z12 = AllocateMemory2D<T>(size/2, size/2);
        T** Z21 = AllocateMemory2D<T>(size/2, size/2);
        T** Z22 = AllocateMemory2D<T>(size/2, size/2);

        T** R1 = AllocateMemory2D<T>(size/2, size/2);
        T** R2 = AllocateMemory2D<T>(size/2, size/2);
        T** R3 = AllocateMemory2D<T>(size/2, size/2);

        split2D<T>(matA, X11, size/2, 0, 0);
        // printMat(X11,size/2,size/2);
        split2D<T>(matA, X12, size/2, 0, size/2);
        split2D<T>(matA, X21, size/2, size/2, 0);
        split2D<T>(matA, X22, size/2, size/2, size/2);
        // printMat(X21,size/2,size/2);

        split2D<T>(matB, Y11, size/2, 0, 0);
        // printMat(Y11,size/2,size/2);
        split2D<T>(matB, Y12, size/2, 0, size/2);
        split2D<T>(matB, Y21, size/2, size/2, 0);
        split2D<T>(matB, Y22, size/2, size/2, size/2);
        // printMat(Y22,size/2,size/2);

        SubMat<T>(X21,X11,R1,size/2);              // T4
        AddMat<T>(Y11,Y12,R2,size/2);              // T9
        matMul_Strassen_v2<T>(R1,R2,R3,size/2);    // Q6
        AddMat<T>(X21,X22,R1,size/2);              // T2
        matMul_Strassen_v2<T>(R1,Y11,R2,size/2);   // Q2
        SubMat<T>(R3,R2,Z22,size/2);               // -U4
        SubMat<T>(Y21,Y11,R1,size/2);              // T8
        matMul_Strassen_v2<T>(X22,R1,R3,size/2);   // Q4
        AddMat<T>(R2,R3,Z21,size/2);               // Z21
        AddMat<T>(X21,X22,R1,size/2);              // T1
        AddMat<T>(Y11,Y22,R2,size/2);              // T6
        matMul_Strassen_v2<T>(R1,R2,Z11,size/2);   // Q1
        AddMat<T>(Z22,Z11,Z22,size/2);             // 
        AddMat<T>(Z11,R3,Z11,size/2);              // U1
        SubMat<T>(Y12,Y22,R1,size/2);              // T7
        matMul_Strassen_v2<T>(X11,R1,R2,size/2);   // Q3
        AddMat<T>(Z22,R2,Z22,size/2);              // Z22
        AddMat<T>(X11,X12,R1,size/2);              // T3
        matMul_Strassen_v2<T>(R1,Y22,R3,size/2);   // Q5
        AddMat<T>(R2,R3,Z12,size/2);               // Z12
        SubMat<T>(Z11,R3,Z11,size/2);              //
        SubMat<T>(X12,X22,R1,size/2);              // T5
        AddMat<T>(Y21,Y22,R2,size/2);              // T10
        matMul_Strassen_v2<T>(R1,R2,R3,size/2);    // Q7
        AddMat<T>(Z11,R3,Z11,size/2);              // Z11
        printMat(Z11,size/2,size/2);

        join2D<T>(Z11,matC,size/2,0,0);
        join2D<T>(Z12,matC,size/2,0,size/2);
        join2D<T>(Z21,matC,size/2,size/2,0);
        join2D<T>(Z22,matC,size/2,size/2,size/2);

        FreeMemory2D<T>(X11);
        FreeMemory2D<T>(X12);
        FreeMemory2D<T>(X21);
        FreeMemory2D<T>(X22);

        FreeMemory2D<T>(Y11);
        FreeMemory2D<T>(Y12);
        FreeMemory2D<T>(Y21);
        FreeMemory2D<T>(Y22);

        FreeMemory2D<T>(Z11);
        FreeMemory2D<T>(Z12);
        FreeMemory2D<T>(Z21);
        FreeMemory2D<T>(Z22);

        FreeMemory2D<T>(R1);
        FreeMemory2D<T>(R2);
        FreeMemory2D<T>(R3);
    }
}