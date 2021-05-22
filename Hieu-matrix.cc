#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>

double* add(double* A, double* B, int size);
double* sub(double* A, double* B, int size);
void split(double* P, double* C, int iB, int jB, int size);
void join(double* C, double* P, int iB, int jB, int size);
double* multiply(double* A, double* B, int n) {
    double* R =  (double*) malloc(sizeof(double)*n*n);
    int n_child = n / 2;
    if (n == 1) {
        R[0] = A[0] * B[0];
    }
    else
    {
        double* A11 = (double*) malloc(sizeof(double)*n_child*n_child);
        double* A12 = (double*) malloc(sizeof(double)*n_child*n_child);
        double* A21 = (double*) malloc(sizeof(double)*n_child*n_child);
        double* A22 = (double*) malloc(sizeof(double)*n_child*n_child);
        double* B11 = (double*) malloc(sizeof(double)*n_child*n_child);
        double* B12 = (double*) malloc(sizeof(double)*n_child*n_child);
        double* B21 = (double*) malloc(sizeof(double)*n_child*n_child);
        double* B22 = (double*) malloc(sizeof(double)*n_child*n_child);

        // Divide matrix A into 4 halves
        split(A, A11, 0, 0, n_child);
        split(A, A12, 0, n / 2, n_child);
        split(A, A21, n / 2, 0, n_child);
        split(A, A22, n / 2, n / 2, n_child);
        /** Dividing matrix B into 4 halves **/
        split(B, B11, 0, 0, n_child);
        split(B, B12, 0, n / 2, n_child);
        split(B, B21, n / 2, 0, n_child);
        split(B, B22, n / 2, n / 2, n_child);

        /**
           * M1 = (A11 + A22)(B11 + B22) M2 = (A21 + A22) B11 M3 = A11 (B12 - B22) M4 =
           * A22 (B21 - B11) M5 = (A11 + A12) B22 M6 = (A21 - A11) (B11 + B12) M7 = (A12 -
           * A22) (B21 + B22)
        **/
        double* M1 = multiply(add(A11, A22, n_child), add(B11, B22, n_child), n_child);
        double* M2 = multiply(add(A21, A22, n_child), B11, n_child);
        double* M3 = multiply(A11, sub(B12, B22, n_child), n_child);
        double* M4 = multiply(A22, sub(B21, B11, n_child), n_child);
        double* M5 = multiply(add(A11, A12, n_child), B22, n_child);
        double* M6 = multiply(sub(A21, A11, n_child), add(B11, B12, n_child), n_child);
        double* M7 = multiply(sub(A12, A22, n_child), add(B21, B22, n_child), n_child);

        /**
          * C11 = M1 + M4 - M5 + M7 C12 = M3 + M5 C21 = M2 + M4 C22 = M1 - M2 + M3 + M6
        **/
        double* C11 = add(sub(add(M1, M4, n_child), M5, n_child), M7, n_child);
        double* C12 = add(M3, M5, n_child);
        double* C21 = add(M2, M4, n_child);
        double* C22 = add(sub(add(M1, M3, n_child), M2, n_child), M6, n_child);

        /** join 4 halves into one result matrix **/
        join(C11, R, 0, 0, n_child);
        join(C12, R, 0, n / 2, n_child);
        join(C21, R, n / 2, 0, n_child);
        join(C22, R, n / 2, n / 2, n_child);
    }
    return R;
}
double* add(double* A, double* B, int size) {
    double* C = (double*) malloc(sizeof(double)*size*size);
    int i,j;
    #pragma omp parallel for default(shared) private(i,j)
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            C[i * size + j] = A[i * size + j] + B[i * size + j];
        }
    }
    return C;
}
double* sub(double* A, double* B, int size) {
    double* C =  (double*) malloc(sizeof(double)*size*size);
    int i,j;
    #pragma omp parallel for default(shared) private(i,j)
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            C[i * size + j] = A[i * size + j] - B[i * size + j];
        }
    }
    return C;
}
void split(double* P, double* C, int iB, int jB, int size)
{
    int i1,i2,j1,j2;
    #pragma omp parallel for default(shared) private(i1,j1)
    for (i1 = 0; i1 < size; i1++) {
        i2=i1+iB;
        for (j1 = 0; j1 < size; j1++) {
            j2=j1+jB;
            C[i1 * size + j1] = P[i2 * size*2 + j2];
        }
    }
}
void join(double* C, double* P, int iB, int jB, int size) {
    int i1,i2,j1,j2;
    #pragma omp parallel for default(shared) private(i1,j1)
    for (i1 = 0; i1 < size; i1++) {
        i2=i1+iB;
        for (j1 = 0; j1 < size; j1++){
            j2=j1+jB;
            P[i2 * size * 2 + j2] = C[i1 * size + j1];
        }
        
    }
}
void print(double* matrix, int size) {
    int i,j;
    #pragma omp parallel for default(shared) private(i,j)
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            printf("%f ", matrix[i * size + j]);
        }
        printf("\n");
    }
    printf("\n");
}
int main(){
    const int n = 256; // This problem scales as n^3. 
                    // This value may need to be adjusted
    double* A = (double*) malloc(sizeof(double)*n*n);
    double* B = (double*) malloc(sizeof(double)*n*n);
    double* C;

    int i, j;
    
    // Cilk Plus array notation
    #pragma omp parallel for default(shared) private(i,j)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i*n+j] = (double)i / (double)n;
            B[i*n+j] = (double)j / (double)n;

        }
    }
    /** Dividing matrix B into 4 halves **/
	double time = omp_get_wtime();
    C = multiply(A, B, n);
  	time = omp_get_wtime() - time;
	
    // popular alg
    // double time = omp_get_wtime();
    // #pragma omp parallel for default(shared) private(i,j,k)  
    // for (i = 0; i < n; i++) {
    //     for (j = 0; j < n; j++) {
    //         for (k = 0; k < n; k++) {
    //             C[i * n + j] += A[i * n + k] * B[k * n + j];
    //         }
    //     }
    // }
    // time = omp_get_wtime() - time;
  
  printf("Checking the results...\n");
  double norm = 0.0;
  #pragma omp parallel for default(shared) private(i,j)
  for (i=0  ; i < n ; i++)
    for (j=0  ; j < n ; j++)
      norm += (C[i*n+j]-(double)(i*j)/(double)n)*(C[i*n+j]-(double)(i*j)/(double)n);
  
  if (norm > 1e-10)
    printf("Something is wrong... Norm is equal to %f\n", norm);
  else
    printf("Yup, we're good!\n");

   printf("Computing time: %f\n", time);  
}
