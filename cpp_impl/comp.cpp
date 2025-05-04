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

template <typename T>
void sparseGEMM_slowed(T *X, int *col_start_pos, int *col_start_neg, int *row_index_pos, int *row_index_neg, T *b, T *Y, int M, int N, int K)
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

/*
 * Called by the driver to register your functions
 * Use add_function(func, description) to add your own functions
 */
void register_functions()
{
    add_function(&sparseGEMM_base, "base_implementation");
    add_function(&sparseGEMM_slowed, "slowed");
}