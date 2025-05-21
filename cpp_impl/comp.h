#ifndef COMP_H
#define COMP_H

#include "common.h"

template <typename T>
void sparseGEMM_csc_base_impl(T *X, const SparseFormatCSC &W_csr, T *b, T *Y, int M, int N, int K)
{
    const int *col_start_pos = W_csr.col_start_pos.data();
    const int *col_start_neg = W_csr.col_start_neg.data();
    const int *row_index_pos = W_csr.row_index_pos.data();
    const int *row_index_neg = W_csr.row_index_neg.data();

    for (int m = 0; m < M; m++)
    {
        for (int n_idx = 0; n_idx < N; n_idx++)
        {
            T y_val = 0;
            for (int k_idx = col_start_pos[n_idx]; k_idx < col_start_pos[n_idx + 1]; k_idx++)
            {
                y_val += X[m * K + row_index_pos[k_idx]];
            }
            for (int k_idx = col_start_neg[n_idx]; k_idx < col_start_neg[n_idx + 1]; k_idx++)
            {
                y_val -= X[m * K + row_index_neg[k_idx]];
            }
            Y[m * N + n_idx] = y_val + b[n_idx];
        }
    }
}

template <typename T, int UNROLL_FACTOR>
void sparseGEMM_csc_unrolled_impl(
    T *X, const SparseFormatCSC &W_csr, T *b, T *Y,
    int M, int N, int K)
{
    const int *col_start_pos = W_csr.col_start_pos.data();
    const int *col_start_neg = W_csr.col_start_neg.data();
    const int *row_index_pos = W_csr.row_index_pos.data();
    const int *row_index_neg = W_csr.row_index_neg.data();

    for (int m = 0; m < M; m++)
    {
        for (int n_idx = 0; n_idx < N; n_idx++)
        {
            T y_pos[UNROLL_FACTOR] = {0};
            T y_neg[UNROLL_FACTOR] = {0};

            int k_pos_loop = col_start_pos[n_idx];
            const int end_pos = col_start_pos[n_idx + 1];

            for (; k_pos_loop + UNROLL_FACTOR <= end_pos; k_pos_loop += UNROLL_FACTOR)
            {
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
                    y_pos[u] += X[m * K + row_index_pos[k_pos_loop + u]];
                }
            }

            T y_pos_final = 0;
            for (int u = 0; u < UNROLL_FACTOR; u++)
            {
                y_pos_final += y_pos[u];
            }

            // remainder loop
            for (; k_pos_loop < end_pos; k_pos_loop++)
            {
                y_pos_final += X[m * K + row_index_pos[k_pos_loop]];
            }

            int k_neg_loop = col_start_neg[n_idx];
            const int end_neg = col_start_neg[n_idx + 1];

            for (; k_neg_loop + UNROLL_FACTOR <= end_neg; k_neg_loop += UNROLL_FACTOR)
            {
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
                    y_neg[u] += X[m * K + row_index_neg[k_neg_loop + u]];
                }
            }

            T y_neg_final = 0;
            for (int u = 0; u < UNROLL_FACTOR; u++)
            {
                y_neg_final += y_neg[u];
            }

            // remainder loop
            for (; k_neg_loop < end_neg; k_neg_loop++)
            {
                y_neg_final += X[m * K + row_index_neg[k_neg_loop]];
            }

            Y[m * N + n_idx] = (y_pos_final - y_neg_final) + b[n_idx];
        }
    }
}

// New sparse GEMM implementation using CSR format for W
template <typename T>
void sparseGEMM_csr_base_impl(T *X, const SparseFormatCSR& W_csr, T *b, T *Y, int M, int N, int K)
{
    // Y is M x N
    // X is M x K
    // W is K x N (num_rows = K, num_cols = N for W_csr)

    // Initialize Y with B values. Y_mn = B_n
    // This needs to be done carefully: Y is M x N, B is N x 1 (or 1 x N, applied to each row of Y)
    for (int m = 0; m < M; ++m) {
        for (int n_col = 0; n_col < N; ++n_col) {
            Y[m * N + n_col] = b[n_col];
        }
    }

    for (int m = 0; m < M; ++m) { // Iterate over rows of X and Y
        for (int k_row = 0; k_row < W_csr.num_rows; ++k_row) { // Iterate over rows of W (which is K)
            T x_val = X[m * K + k_row];

            // Positive contributions from W
            for (int j = W_csr.row_start_pos[k_row]; j < W_csr.row_start_pos[k_row + 1]; ++j) {
                int n_col = W_csr.col_index_pos[j];
                Y[m * N + n_col] += x_val; // W_val is +1
            }

            // Negative contributions from W
            for (int j = W_csr.row_start_neg[k_row]; j < W_csr.row_start_neg[k_row + 1]; ++j) {
                int n_col = W_csr.col_index_neg[j];
                Y[m * N + n_col] -= x_val; // W_val is -1
            }
        }
    }
}

template <typename T, int UNROLL_FACTOR>
void sparseGEMM_csr_unrolled_impl(
    T *X, const SparseFormatCSR& W_csr, T *b, T *Y,
    int M, int N, int K)
{
    // Y is M x N
    // X is M x K
    // W is K x N (num_rows = K, num_cols = N for W_csr)

    // Initialize Y with B values.
    for (int m = 0; m < M; ++m) {
        for (int n_col = 0; n_col < N; ++n_col) {
            Y[m * N + n_col] = b[n_col];
        }
    }

    for (int m = 0; m < M; ++m) { // Iterate over rows of X and Y
        for (int k_row = 0; k_row < W_csr.num_rows; ++k_row) { // Iterate over rows of W (which is K)
            T x_val = X[m * K + k_row];

            // Positive contributions from W
            int j_pos = W_csr.row_start_pos[k_row];
            const int end_pos = W_csr.row_start_pos[k_row + 1];

            // Unrolled loop for positive contributions
            for (; j_pos + UNROLL_FACTOR <= end_pos; j_pos += UNROLL_FACTOR) {
                for (int u = 0; u < UNROLL_FACTOR; ++u) {
                    int n_col = W_csr.col_index_pos[j_pos + u];
                    Y[m * N + n_col] += x_val; // W_val is +1
                }
            }
            // Remainder loop for positive contributions
            for (; j_pos < end_pos; ++j_pos) {
                int n_col = W_csr.col_index_pos[j_pos];
                Y[m * N + n_col] += x_val; // W_val is +1
            }

            // Negative contributions from W
            int j_neg = W_csr.row_start_neg[k_row];
            const int end_neg = W_csr.row_start_neg[k_row + 1];

            // Unrolled loop for negative contributions
            for (; j_neg + UNROLL_FACTOR <= end_neg; j_neg += UNROLL_FACTOR) {
                for (int u = 0; u < UNROLL_FACTOR; ++u) {
                    int n_col = W_csr.col_index_neg[j_neg + u];
                    Y[m * N + n_col] -= x_val; // W_val is -1
                }
            }
            // Remainder loop for negative contributions
            for (; j_neg < end_neg; ++j_neg) {
                int n_col = W_csr.col_index_neg[j_neg];
                Y[m * N + n_col] -= x_val; // W_val is -1
            }
        }
    }
}


#endif // COMP_H

