#ifndef SPARSE_FORMAT_H
#define SPARSE_FORMAT_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

typedef struct ternarySparseFormat_t {
    int *col_start_pos;
    int *col_start_neg;
    int *row_index_pos;
    int *row_index_neg;
    int size;
} ternarySparseFormat;

// Convert ternary matrix to Ternary Sparse Format
ternarySparseFormat *convertTernaryToSparseFormat(int* matrix, int K, int N, int nonZeroPercentage);

bool compare_results(int *result, int *groundTruth, int H, int W);

void GEMM(int *X, int *W, int *b, int *Y, int M, int N, int K);



float *initX(int LEN, int Range);

// Do Sparse GEMM, store results in parameter Y
void sparseGEMM(int* X, ternarySparseFormat *W, int* b, int* Y, int M, int N, int K);

float *generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution);

// Free memory from malloc()
void destroyTernarySparseFormat(ternarySparseFormat *tsf);

#endif // SPARSE_FORMAT_H