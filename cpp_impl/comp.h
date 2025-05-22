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

template <typename T, int TILE_M, int TILE_N_ACCUM, int UNROLL_W_NNZ>
void TCSR_unrolled_tiled(T *X, const TCSRMatrix &W, T *B,
                         T *Y, int M, int N, int K)
{
    const int *row_offsets = W.row_offsets.data();
    const int *encoded_cols = W.encoded_cols.data();
    const int *row_pos_counts = W.row_pos_counts.data();

    T y_tile[TILE_N_ACCUM];

    for (int m_start = 0; m_start < M; m_start += TILE_M)
    {
        int m_end = std::min(m_start + TILE_M, M);

        for (int m = m_start; m < m_end; ++m)
        {
            const T *X_row = X + m * K;
            T *Y_row = Y + m * N;

            for (int n_start = 0; n_start < N; n_start += TILE_N_ACCUM)
            {
                int n_end = std::min(n_start + TILE_N_ACCUM, N);
                int tile_size = n_end - n_start;

                // Initialize tile with bias
                for (int i = 0; i < tile_size; ++i)
                    y_tile[i] = B[n_start + i];

                // Process W rows
                for (int k = 0; k < K; ++k)
                {
                    T x = X_row[k];
                    if (x == static_cast<T>(0))
                        continue;

                    int row_start = row_offsets[k];
                    int pos = row_pos_counts[k];
                    int pos_end = row_start + pos;
                    int row_end = row_offsets[k + 1];

                    // Positive entries
                    for (int i = row_start; i < pos_end;)
                    {
#pragma unroll
                        for (int u = 0; u < UNROLL_W_NNZ && i < pos_end; ++u, ++i)
                        {
                            int n = encoded_cols[i];
                            if (n >= n_start && n < n_end)
                                y_tile[n - n_start] += x;
                        }
                    }

                    // Negative entries
                    for (int i = pos_end; i < row_end;)
                    {
#pragma unroll
                        for (int u = 0; u < UNROLL_W_NNZ && i < row_end; ++u, ++i)
                        {
                            int n = -encoded_cols[i] - 1;
                            if (n >= n_start && n < n_end)
                                y_tile[n - n_start] -= x;
                        }
                    }
                }

                // Store results
                for (int i = 0; i < tile_size; ++i)
                    Y_row[n_start + i] = y_tile[i];
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

// --- Explicit Instantiations ---
// This tells the compiler to generate code for these specific versions in comp.o
template void CSC_base<float>(float *, const SparseFormatCSC &, float *, float *, int, int, int);
template void CSC_base_testing<float>(float *, const SparseFormatCSC &, float *, float *, int, int, int);
template void CCSC_base<float>(float *, const CompressedCSC &, float *, float *, int, int, int);
template void TCSR_base<float>(float *, const TCSRMatrix &, float *, float *, int, int, int);
template void TCSC_base<float>(float *, const TCSCMatrix &, float *, float *, int, int, int);
template void CSC_unrolled<float, 2>(float *, const SparseFormatCSC &, float *, float *, int, int, int);
// If you use other unroll factors or other types for T, you'd add them here.
template void CSC_unrolled<float, 12>(float *, const SparseFormatCSC &, float *, float *, int, int, int);

template void TCSR_unrolled_tiled<float, 12, 8, 8>(float *, const TCSRMatrix &, float *, float *, int, int, int);

template void TCSC_unrolled_tiled<float, 12, 8, 8>(float *, const TCSCMatrix &, float *, float *, int, int, int);

#endif