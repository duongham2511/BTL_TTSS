#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

class Utility {
    public:
    template <typename T>
    static T** AllocateMemory2D(int height, int width)
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
    static void FreeMemory2D(T** arr)
    {
        delete[] *arr;
        delete[] arr;
    }

    static int getPower2(int value)
    {
        int pow = 1;
        while (pow < value) pow = pow*2;
        return pow;
    }

    static void printMat(int** mat, int height, int width)
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

    static void printMat(double** mat, int height, int width)
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

};

#endif