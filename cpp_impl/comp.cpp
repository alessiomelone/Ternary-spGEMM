#include "common.h"

template <typename T>
void sparseGEMM_base(T *X, int *col_start_pos, int *col_start_neg, int *row_index_pos, int *row_index_neg, T *b, T *Y, int M, int N, int K)
{
    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            T y = 0;
            for (int k = col_start_pos[n]; k < col_start_pos[n + 1]; k++)
            {
                y += X[m * K + row_index_pos[k]];
            }
            for (int k = col_start_neg[n]; k < col_start_neg[n + 1]; k++)
            {
                y -= X[m * K + row_index_neg[k]];
            }
            Y[m * N + n] = y + b[n];
        }
    }
}

template <typename T, int UNROLL_FACTOR = 16>
void sparseGEMM_unrolled(
    T *X, int *col_start_pos, int *col_start_neg,
    int *row_index_pos, int *row_index_neg, T *b, T *Y,
    int M, int N, int K)
{

    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            T y_pos[UNROLL_FACTOR] = {0};
            T y_neg[UNROLL_FACTOR] = {0};

            int k_pos = col_start_pos[n];
            const int end_pos = col_start_pos[n + 1];
            
            for (; k_pos + UNROLL_FACTOR <= end_pos; k_pos += UNROLL_FACTOR) {
                #pragma unroll
                for (int u = 0; u < UNROLL_FACTOR; u++) {
                    y_pos[u] += X[m * K + row_index_pos[k_pos + u]];
                }
            }

            T y_pos_final = 0;
            for (int u = 0; u < UNROLL_FACTOR; u++) {
                y_pos_final += y_pos[u];
            }

            // remainder loop
            for (; k_pos < end_pos; k_pos++) {
                y_pos_final += X[m * K + row_index_pos[k_pos]];
            }

            int k_neg = col_start_neg[n];
            const int end_neg = col_start_neg[n + 1];
            
            for (; k_neg + UNROLL_FACTOR <= end_neg; k_neg += UNROLL_FACTOR) {
                #pragma unroll
                for (int u = 0; u < UNROLL_FACTOR; u++) {
                    y_neg[u] += X[m * K + row_index_neg[k_neg + u]];
                }
            }

            T y_neg_final = 0;
            for (int u = 0; u < UNROLL_FACTOR; u++) {
                y_neg_final += y_neg[u];
            }

            // remainder loop
            for (; k_neg < end_neg; k_neg++) {
                y_neg_final += X[m * K + row_index_neg[k_neg]];
            }

            Y[m * N + n] = (y_pos_final - y_neg_final) + b[n];
        }
    }
}

/*
 * Called by the driver to register your functions
 * Use add_function(func, description) to add your own functions
 */
void register_functions()
{
    add_function(&sparseGEMM_base, "base_implementation");
    add_function(&sparseGEMM_unrolled, "unrolled");
}