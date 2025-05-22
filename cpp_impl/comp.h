#ifndef COMP_H
#define COMP_H

#include "common.h"

// Rename and modify sparseGEMM_base to be a specific implementation for SparseFormatCSC
template <typename T>
void CSC_base(T *X, const SparseFormatCSC &W_csr, T *b, T *Y, int M, int N, int K)
{
    const int *col_start_pos = W_csr.col_start_pos.data();
    const int *col_start_neg = W_csr.col_start_neg.data();
    const int *row_index_pos = W_csr.row_index_pos.data();
    const int *row_index_neg = W_csr.row_index_neg.data();

    for (int m = 0; m < M; m++)
    {
        for (int n_idx = 0; n_idx < N; n_idx++) // Renamed n to n_idx
        {
            T y_val = 0;                                                                      // Renamed y to y_val
            for (int k_idx = col_start_pos[n_idx]; k_idx < col_start_pos[n_idx + 1]; k_idx++) // Renamed k to k_idx
            {
                y_val += X[m * K + row_index_pos[k_idx]];
            }
            for (int k_idx = col_start_neg[n_idx]; k_idx < col_start_neg[n_idx + 1]; k_idx++) // Renamed k to k_idx
            {
                y_val -= X[m * K + row_index_neg[k_idx]];
            }
            Y[m * N + n_idx] = y_val + b[n_idx];
        }
    }
}

// comment out X or W memory access to see how much it's slowing it down
template <typename T>
void CSC_base_testing(T *X, const SparseFormatCSC &W_csr, T *b, T *Y, int M, int N, int K)
{
    // grab the column‐pointer arrays once
    const int *col_start_pos = W_csr.col_start_pos.data();
    const int *col_start_neg = W_csr.col_start_neg.data();

    for (int m = 0; m < M; ++m)
    {
        for (int n = 0; n < N; ++n)
        {
            T y_val = 0;

            // simulate the same # of pos‐entries
            int pos_count = col_start_pos[n + 1] - col_start_pos[n];
            for (int i = 0; i < pos_count; ++i)
            {
                y_val += T(1); // replace X[...] and W[...] with constant 1
            }

            // simulate the same # of neg‐entries
            int neg_count = col_start_neg[n + 1] - col_start_neg[n];
            for (int i = 0; i < neg_count; ++i)
            {
                y_val -= T(1); // replace X[...] and W[...] with constant 1
            }

            // keep the bias add so you still write to Y
            Y[m * N + n] = y_val + b[n];
        }
    }
}

template <typename T>
void CCSC_base(T *X, // dense input X, row-major, size M×K
               const CompressedCSC &W,
               T *b,                // bias vector, length N
               T *Y,                // output Y, row-major, size M×N
               int M, int N, int K) // dimensions
{
    // raw pointers into W’s storage
    const int *col_start = W.col_start.data(); // where each column’s blocks begin in vals/row_index
    const int *row_index = W.row_index.data(); // starting row for each encoded block
    const uint8_t *vals = W.vals.data();       // encoded bytes packing 5 ternary values

    for (int m = 0; m < M; ++m) // for each row m of X and Y
    {
        for (int n = 0; n < N; ++n) // for each column n of W and Y
        {
            // five partial sums for the 5 rows inside each block
            T y_val0 = 0;
            T y_val1 = 0;
            T y_val2 = 0;
            T y_val3 = 0;
            T y_val4 = 0;

            // scan through all blocks in column n
            for (int k = col_start[n]; k < col_start[n + 1]; ++k)
            {
                // “row” is the first row index of 5 consecutive rows in W, comprising a block.
                int row = row_index[k];
                const int8_t *d = decodeCCSC[vals[k]]; // decode byte into an array of 5 values (-1/0/1)

                // multiply-add each of the 5 values with X’s corresponding entries
                y_val0 += d[0] * X[m * K + row + 0];
                y_val1 += d[1] * X[m * K + row + 1];
                y_val2 += d[2] * X[m * K + row + 2];
                y_val3 += d[3] * X[m * K + row + 3];
                y_val4 += d[4] * X[m * K + row + 4];
            }

            // combine partials into the full dot product
            T acc = y_val0 + y_val1 + y_val2 + y_val3 + y_val4;

            // add bias for column n and store in Y
            Y[m * N + n] = acc + b[n];
        }
    }
}

template <typename T>
void TCSR_base(T *X, const TCSRMatrix &W, T *B, T *Y, int M, int N, int K)
{
    const int *row_offsets = W.row_offsets.data();
    const int *encoded_cols = W.encoded_cols.data();
    const int *row_pos_counts = W.row_pos_counts.data();

    for (int m = 0; m < M; ++m)
    {
        T *X_row = X + m * K;
        T *Y_row = Y + m * N;

        // Initialize with bias
        for (int n = 0; n < N; ++n)
            Y_row[n] = B[n];

        for (int k = 0; k < K; ++k)
        {
            T x = X_row[k];
            if (x == static_cast<T>(0))
                continue;

            int row_start = row_offsets[k];
            int pos_count = row_pos_counts[k];

            // Positive entries
            for (int i = row_start; i < row_start + pos_count; ++i)
                Y_row[encoded_cols[i]] += x;

            // Negative entries
            for (int i = row_start + pos_count; i < row_offsets[k + 1]; ++i)
                Y_row[-encoded_cols[i] - 1] -= x;
        }
    }
}

template <typename T, int TILE_M, int UNROLL_W_NNZ>
void TCSR_unrolled_tiled(
    T *X_arg,
    const TCSRMatrix &W_tcsr,
    T *B_arg,
    T *Y_arg,
    int M_dim,
    int N_dim,
    int K_dim)
{
    const int *row_offsets_data = W_tcsr.row_offsets.data();
    const int *encoded_cols_data = W_tcsr.encoded_cols.data();
    const int *row_pos_counts_data = W_tcsr.row_pos_counts.data();

    std::vector<T> y_row_accumulator_storage;
    if (N_dim > 0)
    {
        y_row_accumulator_storage.resize(N_dim);
    }
    T *y_row_accumulator = y_row_accumulator_storage.data();

    for (int m_tile_start = 0; m_tile_start < M_dim; m_tile_start += TILE_M)
    {
        const int m_tile_end = std::min(m_tile_start + TILE_M, M_dim);

        for (int m = m_tile_start; m < m_tile_end; ++m)
        {
            const T *X_row_m = X_arg + (size_t)m * K_dim;
            T *Y_row_m = Y_arg + (size_t)m * N_dim;

            for (int n_idx = 0; n_idx < N_dim; ++n_idx)
            {
                y_row_accumulator[n_idx] = B_arg[n_idx];
            }

            for (int k = 0; k < K_dim; ++k)
            {
                const T x_mk_val = X_row_m[k];

                if (x_mk_val == static_cast<T>(0))
                {
                    continue;
                }

                const int W_k_row_start_offset = row_offsets_data[k];
                const int W_k_pos_count = row_pos_counts_data[k];
                const int W_k_pos_entries_end_offset = W_k_row_start_offset + W_k_pos_count;
                const int W_k_row_end_offset = row_offsets_data[k + 1];

                int nnz_ptr_W;

                nnz_ptr_W = W_k_row_start_offset;
                for (; nnz_ptr_W + UNROLL_W_NNZ <= W_k_pos_entries_end_offset; nnz_ptr_W += UNROLL_W_NNZ)
                {
#pragma unroll
                    for (int u = 0; u < UNROLL_W_NNZ; ++u)
                    {
                        y_row_accumulator[encoded_cols_data[nnz_ptr_W + u]] += x_mk_val;
                    }
                }
                for (; nnz_ptr_W < W_k_pos_entries_end_offset; ++nnz_ptr_W)
                {
                    y_row_accumulator[encoded_cols_data[nnz_ptr_W]] += x_mk_val;
                }

                nnz_ptr_W = W_k_pos_entries_end_offset;
                for (; nnz_ptr_W + UNROLL_W_NNZ <= W_k_row_end_offset; nnz_ptr_W += UNROLL_W_NNZ)
                {
#pragma unroll
                    for (int u = 0; u < UNROLL_W_NNZ; ++u)
                    {
                        y_row_accumulator[-encoded_cols_data[nnz_ptr_W + u] - 1] -= x_mk_val;
                    }
                }
                for (; nnz_ptr_W < W_k_row_end_offset; ++nnz_ptr_W)
                {
                    y_row_accumulator[-encoded_cols_data[nnz_ptr_W] - 1] -= x_mk_val;
                }
            }

            for (int n_idx = 0; n_idx < N_dim; ++n_idx)
            {
                Y_row_m[n_idx] = y_row_accumulator[n_idx];
            }
        }
    }
}

template <typename T>
void TCSC_base(T *X, const TCSCMatrix &W_tcsc, T *b, T *Y, int M, int N, int K)
{
    const int *col_offsets = W_tcsc.col_offsets.data();
    const int *encoded_rows = W_tcsc.encoded_rows.data();
    const int *col_pos_counts = W_tcsc.col_pos_counts.data();

    for (int m = 0; m < M; m++)
    {
        T *X_row = X + m * K;
        T *Y_row = Y + m * N;

        for (int n = 0; n < N; n++)
        {
            T y = 0;
            int col_start = col_offsets[n];
            int pos_count = col_pos_counts[n];

            // Positive entries (raw indices)
            for (int k = col_start; k < col_start + pos_count; k++)
            {
                y += X_row[encoded_rows[k]];
            }

            // Negative entries (encoded indices)
            for (int k = col_start + pos_count; k < col_offsets[n + 1]; k++)
            {
                y -= X_row[-encoded_rows[k] - 1];
            }

            Y_row[n] = y + b[n];
        }
    }
}

template <typename T, int TILE_M, int TILE_N, int UNROLL_K>
void TCSC_unrolled_tiled(T *X, const TCSCMatrix &W, T *B,
                         T *Y, int M, int N, int K)
{
    const int *col_offsets = W.col_offsets.data();
    const int *encoded_rows = W.encoded_rows.data();
    const int *col_pos_counts = W.col_pos_counts.data();

    for (int m_start = 0; m_start < M; m_start += TILE_M)
    {
        int m_end = std::min(m_start + TILE_M, M);

        for (int n_start = 0; n_start < N; n_start += TILE_N)
        {
            int n_end = std::min(n_start + TILE_N, N);

            for (int m = m_start; m < m_end; ++m)
            {
                const T *X_row = X + m * K;
                T *Y_row = Y + m * N;

                for (int n = n_start; n < n_end; ++n)
                {
                    T acc = B[n];
                    int off = col_offsets[n];
                    int pos = col_pos_counts[n];
                    int pos_end = off + pos;
                    int col_end = col_offsets[n + 1];

                    // Positive entries
                    int k = off;
                    if (UNROLL_K > 1)
                    {
                        T partials[UNROLL_K] = {0};
#pragma unroll
                        for (; k + UNROLL_K <= pos_end; k += UNROLL_K)
                        {
                            for (int u = 0; u < UNROLL_K; ++u)
                                partials[u] += X_row[encoded_rows[k + u]];
                        }
                        for (int u = 0; u < UNROLL_K; ++u)
                            acc += partials[u];
                    }
                    for (; k < pos_end; ++k)
                        acc += X_row[encoded_rows[k]];

                    // Negative entries
                    k = pos_end;
                    if (UNROLL_K > 1)
                    {
                        T partials[UNROLL_K] = {0};
#pragma unroll
                        for (; k + UNROLL_K <= col_end; k += UNROLL_K)
                        {
                            for (int u = 0; u < UNROLL_K; ++u)
                                partials[u] -= X_row[-encoded_rows[k + u] - 1];
                        }
                        for (int u = 0; u < UNROLL_K; ++u)
                            acc += partials[u];
                    }
                    for (; k < col_end; ++k)
                        acc -= X_row[-encoded_rows[k] - 1];

                    Y_row[n] = acc;
                }
            }
        }
    }
}

// Rename and modify sparseGEMM_unrolled to be a specific implementation for SparseFormatCSC
template <typename T, int UNROLL_FACTOR>
void CSC_unrolled(
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
        { // Renamed n
            T y_pos[UNROLL_FACTOR] = {0};
            T y_neg[UNROLL_FACTOR] = {0};

            int k_pos_loop = col_start_pos[n_idx]; // Renamed k_pos
            const int end_pos = col_start_pos[n_idx + 1];

            for (; k_pos_loop + UNROLL_FACTOR <= end_pos; k_pos_loop += UNROLL_FACTOR)
            {
#pragma unroll
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

            int k_neg_loop = col_start_neg[n_idx]; // Renamed k_neg
            const int end_neg = col_start_neg[n_idx + 1];

            for (; k_neg_loop + UNROLL_FACTOR <= end_neg; k_neg_loop += UNROLL_FACTOR)
            {
#pragma unroll
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

#endif