#ifndef COMP_H
#define COMP_H

#include "common.h"

#ifdef INSTRUMENTATION_RUN
long long flops = 0;
int ds_size = 0;

long long getTotalFlops()
{
    return flops;
}

int getDataStructureSizeInBytes()
{
    return ds_size;
}
#endif

template <typename T>
void BaseCSC(T *X, const BaseTCSC &W_csc, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *col_start_pos = W_csc.col_start_pos.data();
    const int *col_start_neg = W_csc.col_start_neg.data();
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();

    for (int m = 0; m < M; m++)
    {
<<<<<<< HEAD
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
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD:cpp_impl/comp.h
void CCSC_base(T *X, const CompressedCSC &W, T *b, T *Y, int M, int n_col, int N_Rows)
// X: M rows, N_Rows cols
// W: N_Rows rows, n_col cols
// Y: M rows, n_col cols
=======
void CSC_base_testing(T* X, const SparseFormat &W_csr, T *b, T *Y, int M, int N, int K)
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations):cpp_impl/comp.cpp
=======
void CCSC_base(T *X, const CompressedCSC &W, T *b, T *Y, int M, int N, int K)
// X: M rows, K cols
// W: K rows, N cols
// Y: M rows, N cols
>>>>>>> 422cf16 (merge)
=======
void CSC_base_testing(T* X, const SparseFormatCSC &W_csr, T *b, T *Y, int M, int N, int K)
>>>>>>> 6073638 (merge and fix bugs)
=======
void CSC_base_testing(T *X, const SparseFormatCSC &W_csr, T *b, T *Y, int M, int N, int K)
>>>>>>> 9966890 (fix)
=======
void CSC_base_testing(T *X, const BaseTCSC &W_csr, T *b, T *Y, int M, int N, int K)
>>>>>>> 72fe4cb (changed names to BaseTCSC and BaseTCSR)
{
    // grab the column‐pointer arrays once
    const int *col_start_pos = W_csr.col_start_pos.data();
    const int *col_start_neg = W_csr.col_start_neg.data();

    for (int m = 0; m < M; ++m)
    {
        for (int n = 0; n < N; ++n)
=======
        for (int n = 0; n < N; n++)
>>>>>>> 3e9bdeb (added run_benchmark.py)
        {
            T y_val = 0;
            for (int k = col_start_pos[n]; k < col_start_pos[n + 1]; k++)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                y_val += X[m * K + row_index_pos[k]];
            }
            for (int k = col_start_neg[n]; k < col_start_neg[n + 1]; k++)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                y_val -= X[m * K + row_index_neg[k]];
            }
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif
            Y[m * N + n] = y_val + b[n];
        }
    }
}

template <typename T>
void BaseCSC_unroll5(T *X, const BaseTCSC &W_csc,
                     T *b, T *Y,
                     int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *col_start_pos = W_csc.col_start_pos.data();
    const int *col_start_neg = W_csc.col_start_neg.data();
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();

    for (int m = 0; m < M; ++m)
    {
        const T *x_row = X + m * K; // pointer to X[m,0]

        for (int n = 0; n < N; ++n)
        {
            T y_val = 0;

            /* ---------- positive coefficients ---------- */
            int k0 = col_start_pos[n];
            int kEnd = col_start_pos[n + 1];

            /* main body: 5 nz per batch */
            for (; k0 + 4 < kEnd; k0 += 5)
            {
#ifdef INSTRUMENTATION_RUN
                flops += 5;
#endif
                y_val += x_row[row_index_pos[k0]];
                y_val += x_row[row_index_pos[k0 + 1]];
                y_val += x_row[row_index_pos[k0 + 2]];
                y_val += x_row[row_index_pos[k0 + 3]];
                y_val += x_row[row_index_pos[k0 + 4]];
            }
            /* tail */
            for (; k0 < kEnd; ++k0)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                y_val += x_row[row_index_pos[k0]];
            }

            /* ---------- negative coefficients ---------- */
            k0 = col_start_neg[n];
            kEnd = col_start_neg[n + 1];

            for (; k0 + 4 < kEnd; k0 += 5)
            {
                y_val -= x_row[row_index_neg[k0]];
                y_val -= x_row[row_index_neg[k0 + 1]];
                y_val -= x_row[row_index_neg[k0 + 2]];
                y_val -= x_row[row_index_neg[k0 + 3]] y_val -= x_row[row_index_neg[k0 + 4]];
#ifdef INSTRUMENTATION_RUN
                flops += 5;
#endif
            }
            for (; k0 < kEnd; ++k0)
            {
                y_val -= x_row[row_index_neg[k0]];
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }
            /* bias + store */
            Y[m * N + n] = y_val + b[n];
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif
        }
    }
}

template <typename T>
void CCSC_base(T *X, const CompressedCSC &W, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W.getDataStructureSize();
#endif
    const int *col_start = W.col_start.data();
    const int *row_index = W.row_index.data();
    const uint8_t *vals = W.vals.data();

    for (int m = 0; m < M; ++m)
    {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD:cpp_impl/comp.h
        for (int n = 0; n < n_col; n++)
=======
        for (int n = 0; n < N; ++n)        // for each column n of W and Y
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations):cpp_impl/comp.cpp
=======
        for (int n = 0; n < N; n++)
>>>>>>> 422cf16 (merge)
=======
        for (int n = 0; n < N; ++n)        // for each column n of W and Y
>>>>>>> 6073638 (merge and fix bugs)
=======
        for (int n = 0; n < N; ++n) // for each column n of W and Y
>>>>>>> 9966890 (fix)
=======
        for (int n = 0; n < N; ++n)
>>>>>>> 3e9bdeb (added run_benchmark.py)
        {
            T y_val0 = 0;
            T y_val1 = 0;
            T y_val2 = 0;
            T y_val3 = 0;
            T y_val4 = 0;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD:cpp_impl/comp.h
            for (int N_Rows = col_start[n]; N_Rows < col_start[n + 1]; N_Rows++)
=======
            for (int k = col_start[n]; k < col_start[n + 1]; k++)
>>>>>>> 422cf16 (merge)
=======

            // scan through all blocks in column n
            for (int k = col_start[n]; k < col_start[n + 1]; ++k)
>>>>>>> 6073638 (merge and fix bugs)
            {
                // “row” is the first row index of 5 consecutive rows in W, comprising a block.
                int          row = row_index[k];
                const int8_t *d  = decodeCCSC[vals[k]];// decode byte into an array of 5 values (-1/0/1)

                // multiply-add each of the 5 values with X’s corresponding entries
                y_val0 += d[0] * X[m * K + row + 0];
                y_val1 += d[1] * X[m * K + row + 1];
                y_val2 += d[2] * X[m * K + row + 2];
                y_val3 += d[3] * X[m * K + row + 3];
                y_val4 += d[4] * X[m * K + row + 4];
            }

<<<<<<< HEAD
<<<<<<< HEAD
            Y[m * n_col + n + 0] = y_val0 + b[n];
            Y[m * n_col + n + 1] = y_val1 + b[n];
            Y[m * n_col + n + 2] = y_val2 + b[n];
            Y[m * n_col + n + 3] = y_val3 + b[n];
            Y[m * n_col + n + 4] = y_val4 + b[n];
=======

            for (int k = col_start[n]; k < col_start[n + 1]; ++k)
            {
                int row = row_index[k];
                const int8_t *d = decodeCCSC[vals[k]];

                y_val0 += d[0] * X[m * K + row + 0];
                y_val1 += d[1] * X[m * K + row + 1];
                y_val2 += d[2] * X[m * K + row + 2];
                y_val3 += d[3] * X[m * K + row + 3];
                y_val4 += d[4] * X[m * K + row + 4];
#ifdef INSTRUMENTATION_RUN
                flops += 5;
#endif
            }

<<<<<<< HEAD
=======
>>>>>>> 6073638 (merge and fix bugs)
            // combine partials into the full dot product
=======
>>>>>>> 3e9bdeb (added run_benchmark.py)
            T acc = y_val0 + y_val1 + y_val2 + y_val3 + y_val4;
            Y[m * N + n] = acc + b[n];
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations):cpp_impl/comp.cpp
=======
            Y[m * N + n + 0] = y_val0 + b[n];
            Y[m * N + n + 1] = y_val1 + b[n];
            Y[m * N + n + 2] = y_val2 + b[n];
            Y[m * N + n + 3] = y_val3 + b[n];
            Y[m * N + n + 4] = y_val4 + b[n];
>>>>>>> 422cf16 (merge)
=======
>>>>>>> 6073638 (merge and fix bugs)
=======
        #ifdef INSTRUMENTATION_RUN
            flops += 5;
        #endif
>>>>>>> bc77f8b (Add Flop Instrumentation)
=======
#ifdef INSTRUMENTATION_RUN
            flops += 5;
#endif
>>>>>>> da3265c (corrected B in readme)
        }
    }
}

template <typename T>
void TCSR(T *X, const CompressedTCSR &W, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W.getDataStructureSize();
#endif
    const int *row_offsets = W.row_offsets.data();
    const int *encoded_cols = W.encoded_cols.data();
    const int *row_pos_counts = W.row_pos_counts.data();

    for (int m = 0; m < M; ++m)
    {
        T *X_row = X + m * K;
        T *Y_row = Y + m * N;

        for (int n = 0; n < N; ++n)
            Y_row[n] = b[n];

        for (int k = 0; k < K; ++k)
        {
            T x = X_row[k];
            if (x == static_cast<T>(0))
                continue;

            int row_start = row_offsets[k];
            int pos_count = row_pos_counts[k];

            for (int i = row_start; i < row_start + pos_count; ++i)
            {
                Y_row[encoded_cols[i]] += x;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            for (int i = row_start + pos_count; i < row_offsets[k + 1]; ++i)
            {
                Y_row[-encoded_cols[i] - 1] -= x;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }
        }
    }
}

template <typename T, int TILE_M, int UNROLL_W_NNZ>
void TCSR_unrolled_tiled(T *X, const CompressedTCSR &W, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W.getDataStructureSize();
#endif
    const int *row_offsets = W.row_offsets.data();
    const int *encoded_cols = W.encoded_cols.data();
    const int *row_pos_counts = W.row_pos_counts.data();

    std::vector<T> y_row_accumulator_storage;
    if (N > 0)
    {
        y_row_accumulator_storage.resize(N); // TODO: Why is this not in the init code ?
    }
    T *y_row_accumulator = y_row_accumulator_storage.data();

    for (int m_tile_start = 0; m_tile_start < M; m_tile_start += TILE_M)
    {
        const int m_tile_end = std::min(m_tile_start + TILE_M, M);

        for (int m = m_tile_start; m < m_tile_end; ++m)
        {
            const T *X_row_m = X + (size_t)m * K;
            T *Y_row_m = Y + (size_t)m * N;

            for (int n = 0; n < N; ++n)
            {
                y_row_accumulator[n] = b[n];
            }

            for (int k = 0; k < K; ++k)
            {
                const T x_mk_val = X_row_m[k];

                if (x_mk_val == static_cast<T>(0))
                {
                    continue;
                }

                const int W_k_row_start_offset = row_offsets[k];
                const int W_k_pos_count = row_pos_counts[k];
                const int W_k_pos_entries_end_offset = W_k_row_start_offset + W_k_pos_count;
                const int W_k_row_end_offset = row_offsets[k + 1];

                int nnz_ptr_W;

                nnz_ptr_W = W_k_row_start_offset;
                for (; nnz_ptr_W + UNROLL_W_NNZ <= W_k_pos_entries_end_offset; nnz_ptr_W += UNROLL_W_NNZ)
                {
                    for (int u = 0; u < UNROLL_W_NNZ; ++u)
                    {
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                        y_row_accumulator[encoded_cols[nnz_ptr_W + u]] += x_mk_val;
                    }
                }
                for (; nnz_ptr_W < W_k_pos_entries_end_offset; ++nnz_ptr_W)
                {
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                    y_row_accumulator[encoded_cols[nnz_ptr_W]] += x_mk_val;
                }

                nnz_ptr_W = W_k_pos_entries_end_offset;
                for (; nnz_ptr_W + UNROLL_W_NNZ <= W_k_row_end_offset; nnz_ptr_W += UNROLL_W_NNZ)
                {
                    for (int u = 0; u < UNROLL_W_NNZ; ++u)
                    {
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                        y_row_accumulator[-encoded_cols[nnz_ptr_W + u] - 1] -= x_mk_val;
                    }
                }
                for (; nnz_ptr_W < W_k_row_end_offset; ++nnz_ptr_W)
                {
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                    y_row_accumulator[-encoded_cols[nnz_ptr_W] - 1] -= x_mk_val;
                }
            }

            for (int n = 0; n < N; ++n)
            {
                Y_row_m[n] = y_row_accumulator[n];
            }
        }
    }
}

template <typename T>
void TCSC(T *X, const CompressedTCSC &W, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W.getDataStructureSize();
#endif
    const int *col_offsets = W.col_offsets.data();
    const int *encoded_rows = W.encoded_rows.data();
    const int *col_pos_counts = W.col_pos_counts.data();

    for (int m = 0; m < M; m++)
    {
        T *X_row = X + m * K;
        T *Y_row = Y + m * N;

        for (int n = 0; n < N; n++)
        {
            T y = 0;
            int col_start = col_offsets[n];
            int pos_count = col_pos_counts[n];

            for (int k = col_start; k < col_start + pos_count; k++)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                y += X_row[encoded_rows[k]];
            }

            for (int k = col_start + pos_count; k < col_offsets[n + 1]; k++)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                y -= X_row[-encoded_rows[k] - 1];
            }
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif
            Y_row[n] = y + b[n];
        }
    }
}

template <typename T, int TILE_M, int TILE_N, int UNROLL_K>
void TCSC_unrolled_tiled(T *X, const CompressedTCSC &W, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W.getDataStructureSize();
#endif
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
                    T acc = b[n];
                    int off = col_offsets[n];
                    int pos = col_pos_counts[n];
                    int pos_end = off + pos;
                    int col_end = col_offsets[n + 1];

                    int k = off;
                    if (UNROLL_K > 1)
                    {
                        T partials[UNROLL_K] = {0};
                        for (; k + UNROLL_K <= pos_end; k += UNROLL_K)
                        {
                            for (int u = 0; u < UNROLL_K; ++u)
                            {
                                partials[u] += X_row[encoded_rows[k + u]];
#ifdef INSTRUMENTATION_RUN
                                flops++;
#endif
                            }
                        }
                        for (int u = 0; u < UNROLL_K; ++u)
                        {
                            acc += partials[u];
#ifdef INSTRUMENTATION_RUN
                            flops++;
#endif
                        }
                    }
                    for (; k < pos_end; ++k)
                    {
                        acc += X_row[encoded_rows[k]];
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                    }

                    k = pos_end;
                    if (UNROLL_K > 1)
                    {
                        T partials[UNROLL_K] = {0};
                        for (; k + UNROLL_K <= col_end; k += UNROLL_K)
                        {
                            for (int u = 0; u < UNROLL_K; ++u)
                            {
                                partials[u] -= X_row[-encoded_rows[k + u] - 1];
#ifdef INSTRUMENTATION_RUN
                                flops++;
#endif
                            }
                        }
                        for (int u = 0; u < UNROLL_K; ++u)
                        {
                            acc += partials[u];
#ifdef INSTRUMENTATION_RUN
                            flops++;
#endif
                        }
                    }
                    for (; k < col_end; ++k)
                    {
                        acc -= X_row[-encoded_rows[k] - 1];
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                    }
                    Y_row[n] = acc;
                }
            }
        }
    }
}

template <typename T, int UNROLL_FACTOR>
void BaseCSC_unr(T *X, const BaseTCSC &W_csc, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *col_start_pos = W_csc.col_start_pos.data();
    const int *col_start_neg = W_csc.col_start_neg.data();
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();

    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            T y_pos[UNROLL_FACTOR] = {0};
            T y_neg[UNROLL_FACTOR] = {0};

            int k_pos_loop = col_start_pos[n];
            const int end_pos = col_start_pos[n + 1];

            for (; k_pos_loop + UNROLL_FACTOR <= end_pos; k_pos_loop += UNROLL_FACTOR)
            {
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
                    y_pos[u] += X[m * K + row_index_pos[k_pos_loop + u]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
            }

            T y_pos_final = 0;
            for (int u = 0; u < UNROLL_FACTOR; u++)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                y_pos_final += y_pos[u];
            }

            for (; k_pos_loop < end_pos; k_pos_loop++)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                y_pos_final += X[m * K + row_index_pos[k_pos_loop]];
            }

            int k_neg_loop = col_start_neg[n];
            const int end_neg = col_start_neg[n + 1];

            for (; k_neg_loop + UNROLL_FACTOR <= end_neg; k_neg_loop += UNROLL_FACTOR)
            {
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                    y_neg[u] += X[m * K + row_index_neg[k_neg_loop + u]];
                }
            }

            T y_neg_final = 0;
            for (int u = 0; u < UNROLL_FACTOR; u++)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                y_neg_final += y_neg[u];
            }

            for (; k_neg_loop < end_neg; k_neg_loop++)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                y_neg_final += X[m * K + row_index_neg[k_neg_loop]];
            }
#ifdef INSTRUMENTATION_RUN
            flops += 2;
#endif
            Y[m * N + n] = (y_pos_final - y_neg_final) + b[n];
        }
    }
}

template <typename T>
void BaseCSR(T *X, const BaseTCSR &W_csr, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csr.getDataStructureSize();
#endif
    const int *row_start_pos = W_csr.row_start_pos.data();
    const int *row_start_neg = W_csr.row_start_neg.data();
    const int *col_index_pos = W_csr.col_index_pos.data();
    const int *col_index_neg = W_csr.col_index_neg.data();

    for (int m = 0; m < M; ++m)
    {
        for (int n = 0; n < N; ++n)
        {
            Y[m * N + n] = b[n];
        }
    }

    for (int m = 0; m < M; ++m)
    {
        for (int k = 0; k < K; ++k)
        {
            T x_val = X[m * K + k];
            if (x_val == static_cast<T>(0))
                continue;

            for (int j = row_start_pos[k]; j < row_start_pos[k + 1]; ++j)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                int n_col = col_index_pos[j];
                Y[m * N + n_col] += x_val;
            }

            for (int j = row_start_neg[k]; j < row_start_neg[k + 1]; ++j)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                int n_col = col_index_neg[j];
                Y[m * N + n_col] -= x_val;
            }
        }
    }
}

template <typename T, int UNROLL_FACTOR>
void BaseCSR_unr(T *X, const BaseTCSR &W_csr, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csr.getDataStructureSize();
#endif
    const int *row_start_pos = W_csr.row_start_pos.data();
    const int *row_start_neg = W_csr.row_start_neg.data();
    const int *col_index_pos = W_csr.col_index_pos.data();
    const int *col_index_neg = W_csr.col_index_neg.data();

    for (int m = 0; m < M; ++m)
    {
        for (int n = 0; n < N; ++n)
        {
            Y[m * N + n] = b[n];
        }
    }

    for (int m = 0; m < M; ++m)
    {
        for (int k = 0; k < K; ++k)
        {
            T x_val = X[m * K + k];
            if (x_val == static_cast<T>(0))
                continue;

            int j_pos = row_start_pos[k];
            const int end_pos = row_start_pos[k + 1];

            for (; j_pos + UNROLL_FACTOR <= end_pos; j_pos += UNROLL_FACTOR)
            {
                for (int u = 0; u < UNROLL_FACTOR; ++u)
                {
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                    int n_col = col_index_pos[j_pos + u];
                    Y[m * N + n_col] += x_val;
                }
            }
            for (; j_pos < end_pos; ++j_pos)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                int n_col = col_index_pos[j_pos];
                Y[m * N + n_col] += x_val;
            }

            int j_neg = row_start_neg[k];
            const int end_neg = row_start_neg[k + 1];

            for (; j_neg + UNROLL_FACTOR <= end_neg; j_neg += UNROLL_FACTOR)
            {
                for (int u = 0; u < UNROLL_FACTOR; ++u)
                {
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                    int n_col = col_index_neg[j_neg + u];
                    Y[m * N + n_col] -= x_val;
                }
            }
            for (; j_neg < end_neg; ++j_neg)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                int n_col = col_index_neg[j_neg];
                Y[m * N + n_col] -= x_val;
            }
        }
    }
}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD:cpp_impl/comp.h
// Base implementations
template <typename T>
void sparseGEMM_csc_base_impl(T *X, const SparseFormatCSC &W_csc, T *b, T *Y, int M, int n_col, int N_Rows)
=======
template <typename T, int TILE_M, int TILE_N, int UNROLL_FACTOR>
void BaseCSC_unr_tiled(T *X, const BaseTCSC &W_csc, T *b, T *Y, int M, int N, int K)
>>>>>>> 5ab0bc5 (added 1000x4096x16384 and 4000x4096x16384 to run_benchmark.py)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *col_start_pos = W_csc.col_start_pos.data();
    const int *col_start_neg = W_csc.col_start_neg.data();
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();
<<<<<<< HEAD
=======
// --- Explicit Instantiations ---
// This tells the compiler to generate code for these specific versions in comp.o
<<<<<<< HEAD
template void CSC_base<float>(float *, const SparseFormat &, float *, float *, int, int, int);
template void CSC_base_testing<float>(float *, const SparseFormat &, float *, float *, int, int, int);
=======
// --- Explicit Instantiations ---
// This tells the compiler to generate code for these specific versions in comp.o
template void CSC_base<float>(float *, const SparseFormat &, float *, float *, int, int, int);
>>>>>>> 422cf16 (merge)
=======
template void CSC_base<float>(float *, const SparseFormatCSC &, float *, float *, int, int, int);
template void CSC_base_testing<float>(float *, const SparseFormatCSC &, float *, float *, int, int, int);
>>>>>>> 6073638 (merge and fix bugs)
=======
// --- Explicit Instantiations ---
// This tells the compiler to generate code for these specific versions in comp.o
template void CSC_base<float>(float *, const SparseFormatCSC &, float *, float *, int, int, int);
template void CSC_base_testing<float>(float *, const SparseFormatCSC &, float *, float *, int, int, int);
>>>>>>> f2770dd (TCSC tiled is 1.5)
template void CCSC_base<float>(float *, const CompressedCSC &, float *, float *, int, int, int);
template void TCSR_base<float>(float *, const TCSRMatrix &, float *, float *, int, int, int);
template void TCSC_base<float>(float *, const TCSCMatrix &, float *, float *, int, int, int);
template void CSC_unrolled<float, 2>(float *, const SparseFormatCSC &, float *, float *, int, int, int);
// If you use other unroll factors or other types for T, you'd add them here.
<<<<<<< HEAD
<<<<<<< HEAD
template void CSC_unrolled<float, 12>(float *, const SparseFormat &, float *, float *, int, int, int);
<<<<<<< HEAD
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations):cpp_impl/comp.cpp
=======
>>>>>>> 422cf16 (merge)
=======
template void CSC_unrolled<float, 12>(float *, const SparseFormatCSC &, float *, float *, int, int, int);
>>>>>>> 6073638 (merge and fix bugs)

template void TCSR_unrolled<float, 12>(float *, const TCSRMatrix &, float *, float *, int, int, int);
template void TCSC_unrolled<float, 12>(float *, const TCSCMatrix &, float *, float *, int, int, int);

template void TCSC_unrolled_tiled<float, 12, 32, 32>(float *, const TCSCMatrix &, float *, float *, int, int, int);

=======
>>>>>>> 9966890 (fix)
=======
template void CSC_unrolled<float, 12>(float *, const SparseFormatCSC &, float *, float *, int, int, int);

template void TCSR_unrolled_tiled<float, 12, 8, 8>(float *, const TCSRMatrix &, float *, float *, int, int, int);

template void TCSC_unrolled_tiled<float, 12, 8, 8>(float *, const TCSCMatrix &, float *, float *, int, int, int);

>>>>>>> f2770dd (TCSC tiled is 1.5)
=======
>>>>>>> 84e2726 (Tiling for TCSR improved but still needs works)
=======

    for (int m_start = 0; m_start < M; m_start += TILE_M)
    {
        int m_end = std::min(m_start + TILE_M, M);

        for (int n_start = 0; n_start < N; n_start += TILE_N)
        {
            int n_end = std::min(n_start + TILE_N, N);

            for (int m = m_start; m < m_end; m++)
            {
                for (int n = n_start; n < n_end; n++)
                {
                    T y_pos[UNROLL_FACTOR] = {0};
                    T y_neg[UNROLL_FACTOR] = {0};

                    int k_pos_loop = col_start_pos[n];
                    const int end_pos = col_start_pos[n + 1];

                    for (; k_pos_loop + UNROLL_FACTOR <= end_pos; k_pos_loop += UNROLL_FACTOR)
                    {
                        for (int u = 0; u < UNROLL_FACTOR; u++)
                        {
#ifdef INSTRUMENTATION_RUN
                            flops++;
#endif
                            y_pos[u] += X[m * K + row_index_pos[k_pos_loop + u]];
                        }
                    }

                    T y_pos_final = 0;
                    for (int u = 0; u < UNROLL_FACTOR; u++)
                    {
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                        y_pos_final += y_pos[u];
                    }

                    for (; k_pos_loop < end_pos; k_pos_loop++)
                    {
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                        y_pos_final += X[m * K + row_index_pos[k_pos_loop]];
                    }

                    int k_neg_loop = col_start_neg[n];
                    const int end_neg = col_start_neg[n + 1];

                    for (; k_neg_loop + UNROLL_FACTOR <= end_neg; k_neg_loop += UNROLL_FACTOR)
                    {
                        for (int u = 0; u < UNROLL_FACTOR; u++)
                        {
#ifdef INSTRUMENTATION_RUN
                            flops++;
#endif
                            y_neg[u] += X[m * K + row_index_neg[k_neg_loop + u]];
                        }
                    }

                    T y_neg_final = 0;
                    for (int u = 0; u < UNROLL_FACTOR; u++)
                    {
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                        y_neg_final += y_neg[u];
                    }

                    for (; k_neg_loop < end_neg; k_neg_loop++)
                    {
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                        y_neg_final += X[m * K + row_index_neg[k_neg_loop]];
                    }
#ifdef INSTRUMENTATION_RUN
                    flops += 2;
#endif
                    Y[m * N + n] = (y_pos_final - y_neg_final) + b[n];
                }
            }
        }
    }
}

>>>>>>> 5ab0bc5 (added 1000x4096x16384 and 4000x4096x16384 to run_benchmark.py)
#endif