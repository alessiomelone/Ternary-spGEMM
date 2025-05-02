#include "sparse_format.h"

// Do Sparse GEMM, store results in parameter Y
void sparseGEMM(int* X, ternarySparseFormat *W, int* b, int* Y, int M, int N, int K) {
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            int y = 0;
            for (int k = W->col_start_pos[n]; k < W->col_start_pos[n + 1]; k++) {
                y += X[m * K + W->row_index_pos[k]];
            }
            for (int k = W->col_start_neg[n]; k < W->col_start_neg[n + 1]; k++) {
                y -= X[m * K + W->row_index_neg[k]];
            }
            Y[m * N + n] = y + b[n];
        }
    }
}

// Convert ternary matrix to Ternary Sparse Format
ternarySparseFormat *convertTernaryToSparseFormat(int* matrix, int K, int N, int nonZeroPercentage) {
    // TODO: Verify sizes of each sub array
    int nonZeroVals = (K * N) / (double) (nonZeroPercentage) + 1;
    ternarySparseFormat *tsf = malloc(sizeof(ternarySparseFormat));
    tsf->size = N+1;
    tsf->col_start_pos = malloc(tsf->size * sizeof(int));
    tsf->col_start_neg = malloc(tsf->size * sizeof(int));
    tsf->row_index_pos = malloc(nonZeroVals * sizeof(int));
    tsf->row_index_neg = malloc(nonZeroVals * sizeof(int));

    int column_start_pos = 0;
    int column_start_neg = 0;
    int row_index_pos_ind = 0;
    int row_index_neg_ind = 0;
    int n;
    for (n = 0; n < N; n++) {
        tsf->col_start_pos[n] = column_start_pos;
        tsf->col_start_neg[n] = column_start_neg;
        for (int k = 0; k < K; k++) {
            if (matrix[k * N + n] >= 1) {
                column_start_pos++;
                tsf->row_index_pos[row_index_pos_ind++] = k;
            }
            else if (matrix[k * N + n] <= -1) {
                column_start_neg++;
                tsf->row_index_neg[row_index_neg_ind++] = k;
            }
        }
    }
    tsf->col_start_pos[n] = column_start_pos;
    tsf->col_start_neg[n] = column_start_neg;
    return tsf;
}

// Free memory from malloc()
void destroyTernarySparceFormat(ternarySparseFormat *tsf) {
    free(tsf->col_start_pos);
    free(tsf->col_start_neg);
    free(tsf->row_index_pos);
    free(tsf->row_index_neg);
    free(tsf);
}