#ifndef SPARSE_FORMAT_H
#define SPARSE_FORMAT_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#define ERROR_TOLERANCE ((float) (10e-2))

typedef struct ternarySparseFormat_t {
    int *col_start_pos;
    int *col_start_neg;
    int *row_index_pos;
    int *row_index_neg;
    int size;
} ternarySparseFormat;

// Convert ternary matrix to Ternary Sparse Format
ternarySparseFormat *convertTernaryToSparseFormat(int* matrix, int K, int N, int nonZeroPercentage);

bool compare_results(float *result, float *groundTruth, int H, int W);

void GEMM(float *X, int *W, float *b, float *Y, int M, int N, int K);

float *initX(int LEN, int Range);

// Do Sparse GEMM, store results in parameter Y
void sparseGEMM(float* X, ternarySparseFormat *W, float* b, float* Y, int M, int N, int K);

int *generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution);

// Free memory from malloc()
void destroyTernarySparseFormat(ternarySparseFormat *tsf);

#endif // SPARSE_FORMAT_H