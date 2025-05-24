#ifndef COMP_H
#define COMP_H

#include "common.h"
#include "data_structures/BlockedTCSC.h"
#include <iostream>

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
        for (int n = 0; n < N; n++)
        {
            T y_val = 0;

            // Process positive values
            for (int k = col_start_pos[n]; k < col_start_pos[n + 1]; k++)
            {
                T x_val = X[m * K + row_index_pos[k]];
                y_val += x_val;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            // Process negative values
            for (int k = col_start_neg[n]; k < col_start_neg[n + 1]; k++)
            {
                T x_val = X[m * K + row_index_neg[k]];
                y_val -= x_val;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

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
    const int *col_start_5 = W.col_start_5.data();
    const int *row_index_5 = W.row_index_5.data();

    const int *row_index_pos = W.row_index_pos.data();
    const int *row_index_neg = W.row_index_neg.data();

    const int *col_start_pos = W.col_start_pos.data();
    const int *col_start_neg = W.col_start_neg.data();

    const uint8_t *vals_5 = W.vals_5.data();

    for (int m = 0; m < M; ++m)
    {
        for (int n = 0; n < N; ++n)
        {
            // printf("n=%d, m=%d\n", n, m);
            T y_val0 = 0;
            T y_val1 = 0;
            T y_val2 = 0;
            T y_val3 = 0;
            T y_val4 = 0;

            // blocks of 5
            for (int k = col_start_5[n]; k < col_start_5[n + 1]; ++k)
            {
                int row = row_index_5[k];
                const int8_t *d = decode5[vals_5[k]];

                y_val0 += d[0] * X[m * K + row + 0];
                y_val1 += d[1] * X[m * K + row + 1];
                y_val2 += d[2] * X[m * K + row + 2];
                y_val3 += d[3] * X[m * K + row + 3];
                y_val4 += d[4] * X[m * K + row + 4];
#ifdef INSTRUMENTATION_RUN
                flops += 5;
#endif
            }

            // positive
            for (int k = col_start_pos[n]; k < col_start_pos[n + 1]; ++k)
            {
                int row = row_index_pos[k];
                y_val0 += X[m * K + row];
            }
#ifdef INSTRUMENTATION_RUN
            flops += col_start_pos[n + 1] - col_start_pos[n];
#endif

            // negative
            for (int k = col_start_neg[n]; k < col_start_neg[n + 1]; ++k)
            {
                int row = row_index_neg[k];
                y_val0 -= X[m * K + row];
            }
#ifdef INSTRUMENTATION_RUN
            flops += col_start_neg[n + 1] - col_start_neg[n];
#endif

            T acc = y_val0 + y_val1 + y_val2 + y_val3 + y_val4;
            Y[m * N + n] = acc + b[n];
#ifdef INSTRUMENTATION_RUN
            flops += 5;
#endif
        }
    }
}

template <typename T>
void CCSC_unr(T *X, const CompressedCSC &W, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W.getDataStructureSize();
#endif
    const int *col_start_5 = W.col_start_5.data();
    const int *row_index_5 = W.row_index_5.data();

    const int *row_index_pos = W.row_index_pos.data();
    const int *row_index_neg = W.row_index_neg.data();

    const int *col_start_pos = W.col_start_pos.data();
    const int *col_start_neg = W.col_start_neg.data();

    const uint8_t *vals_5 = W.vals_5.data();

    for (int m = 0; m < M; ++m)
    {
        for (int n = 0; n < N; ++n)
        {
            // printf("n=%d, m=%d\n", n, m);
            T y_val0 = 0;
            T y_val1 = 0;
            T y_val2 = 0;
            T y_val3 = 0;
            T y_val4 = 0;

            // blocks of 5
            for (int k = col_start_5[n]; k < col_start_5[n + 1]; ++k)
            {
                int row = row_index_5[k];
                const int8_t *d = decode5[vals_5[k]];

                y_val0 += d[0] * X[m * K + row + 0];
                y_val1 += d[1] * X[m * K + row + 1];
                y_val2 += d[2] * X[m * K + row + 2];
                y_val3 += d[3] * X[m * K + row + 3];
                y_val4 += d[4] * X[m * K + row + 4];
#ifdef INSTRUMENTATION_RUN
                flops += 5;
#endif
            }

            // positive (fully unrolled by 5)
            T y_pos0 = 0, y_pos1 = 0, y_pos2 = 0, y_pos3 = 0, y_pos4 = 0;

            int k_pos_loop = col_start_pos[n];
            const int end_pos = col_start_pos[n + 1];

            for (; k_pos_loop + 5 <= end_pos; k_pos_loop += 5)
            {
                int row0 = row_index_pos[k_pos_loop];
                int row1 = row_index_pos[k_pos_loop + 1];
                int row2 = row_index_pos[k_pos_loop + 2];
                int row3 = row_index_pos[k_pos_loop + 3];
                int row4 = row_index_pos[k_pos_loop + 4];

                y_pos0 += X[m * K + row0];
                y_pos1 += X[m * K + row1];
                y_pos2 += X[m * K + row2];
                y_pos3 += X[m * K + row3];
                y_pos4 += X[m * K + row4];
#ifdef INSTRUMENTATION_RUN
                flops += 5;
#endif
            }

            T y_pos_final = y_pos0 + y_pos1 + y_pos2 + y_pos3 + y_pos4;
#ifdef INSTRUMENTATION_RUN
            flops += 4; // reductions
#endif

            for (; k_pos_loop < end_pos; ++k_pos_loop)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                int row = row_index_pos[k_pos_loop];
                y_pos_final += X[m * K + row];
            }

            y_val0 += y_pos_final;

            // negative (fully unrolled by 5)
            T y_neg0 = 0, y_neg1 = 0, y_neg2 = 0, y_neg3 = 0, y_neg4 = 0;

            int k_neg_loop = col_start_neg[n];
            const int end_neg = col_start_neg[n + 1];

            for (; k_neg_loop + 5 <= end_neg; k_neg_loop += 5)
            {
                int row0 = row_index_neg[k_neg_loop];
                int row1 = row_index_neg[k_neg_loop + 1];
                int row2 = row_index_neg[k_neg_loop + 2];
                int row3 = row_index_neg[k_neg_loop + 3];
                int row4 = row_index_neg[k_neg_loop + 4];

                y_neg0 += X[m * K + row0];
                y_neg1 += X[m * K + row1];
                y_neg2 += X[m * K + row2];
                y_neg3 += X[m * K + row3];
                y_neg4 += X[m * K + row4];
#ifdef INSTRUMENTATION_RUN
                flops += 5;
#endif
            }

            T y_neg_final = y_neg0 + y_neg1 + y_neg2 + y_neg3 + y_neg4;
#ifdef INSTRUMENTATION_RUN
            flops += 4; // reductions
#endif

            for (; k_neg_loop < end_neg; ++k_neg_loop)
            {
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
                int row = row_index_neg[k_neg_loop];
                y_neg_final += X[m * K + row];
            }

            y_val0 -= y_neg_final;

            T acc = y_val0 + y_val1 + y_val2 + y_val3 + y_val4;
            Y[m * N + n] = acc + b[n];
#ifdef INSTRUMENTATION_RUN
            flops += 5;
#endif
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

            // Process positive values
            for (int k = col_start; k < col_start + pos_count; k++)
            {
                T x_val = X_row[encoded_rows[k]];
                y += x_val;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            // Process negative values
            for (int k = col_start + pos_count; k < col_offsets[n + 1]; k++)
            {
                T x_val = X_row[-encoded_rows[k] - 1];
                y -= x_val;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            Y_row[n] = y + b[n];
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif
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

    // Initialize Y with bias
    for (int m = 0; m < M; ++m)
    {
        for (int n = 0; n < N; ++n)
        {
            Y[m * N + n] = b[n];
        }
    }

    // For each row in X
    for (int m = 0; m < M; ++m)
    {
        T *X_row = X + m * K;
        T *Y_row = Y + m * N;

        // For each row in W
        for (int k = 0; k < K; ++k)
        {
            T x_val = X_row[k];
            if (x_val != static_cast<T>(0))
            {
                // Process positive values
                for (int j = row_start_pos[k]; j < row_start_pos[k + 1]; ++j)
                {
                    Y_row[col_index_pos[j]] += x_val;
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                // Process negative values
                for (int j = row_start_neg[k]; j < row_start_neg[k + 1]; ++j)
                {
                    Y_row[col_index_neg[j]] -= x_val;
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
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

template <typename T, int TILE_M, int TILE_N, int UNROLL_FACTOR>
void BaseCSC_unr_tiled(T *X, const BaseTCSC &W_csc, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *col_start_pos = W_csc.col_start_pos.data();
    const int *col_start_neg = W_csc.col_start_neg.data();
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();

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

template <typename T, int B>
void BlockedCSC(T *X, const BlockedTCSC<B> &W_csc, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *col_start_pos = W_csc.col_start_pos.data();
    const int *col_start_neg = W_csc.col_start_neg.data();
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();

    // Process each row of Y
    for (int m = 0; m < M; m++)
    {
        // Process each block of K
        for (int k_block = 0; k_block < K / B; k_block++)
        {
            for (int n = k_block * N; n < k_block * N + N; n++)
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
                Y[m * N + n % N] += y;
            }
        }

        // Add bias after processing the entire row
        for (int n = 0; n < N; n++)
        {
            Y[m * N + n] += b[n];
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif
        }
    }
}

template <typename T, int B>
void BlockedCSC_unr4(T *X, const BlockedTCSC<B> &W_csc, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *col_start_pos = W_csc.col_start_pos.data();
    const int *col_start_neg = W_csc.col_start_neg.data();
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();

    // Process each row of Y
    for (int m = 0; m < M; m++)
    {
        // Process each block of K
        for (int k_block = 0; k_block < K / B; k_block++)
        {
            for (int n = k_block * N; n < k_block * N + N; n++)
            {
                T y_pos[4] = {0};
                T y_neg[4] = {0};
                int k_pos = col_start_pos[n];
                int k_neg = col_start_neg[n];
                const int end_pos = col_start_pos[n + 1];
                const int end_neg = col_start_neg[n + 1];

                // Process positive values with unrolling
                for (; k_pos + 4 <= end_pos; k_pos += 4)
                {
                    for (int u = 0; u < 4; u++)
                    {
                        y_pos[u] += X[m * K + row_index_pos[k_pos + u]];
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                    }
                }
                // Handle remaining positive values
                T y_pos_final = 0;
                for (int u = 0; u < 4; u++)
                {
                    y_pos_final += y_pos[u];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
                for (; k_pos < end_pos; k_pos++)
                {
                    y_pos_final += X[m * K + row_index_pos[k_pos]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                // Process negative values with unrolling
                for (; k_neg + 4 <= end_neg; k_neg += 4)
                {
                    for (int u = 0; u < 4; u++)
                    {
                        y_neg[u] += X[m * K + row_index_neg[k_neg + u]];
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                    }
                }
                // Handle remaining negative values
                T y_neg_final = 0;
                for (int u = 0; u < 4; u++)
                {
                    y_neg_final += y_neg[u];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
                for (; k_neg < end_neg; k_neg++)
                {
                    y_neg_final += X[m * K + row_index_neg[k_neg]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                Y[m * N + n % N] += (y_pos_final - y_neg_final);
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }
        }

        // Add bias after processing the entire row
        for (int n = 0; n < N; n++)
        {
            Y[m * N + n] += b[n];
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif
        }
    }
}

#endif