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
void ICSR_base(const T *X, const ICSR &W_icsr, const T *B, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_icsr.getDataStructureSize();
#endif

    const int *const row_offsets_data = W_icsr.row_offsets.data();
    const int *const encoded_cols_data = W_icsr.encoded_cols.data();

    for (int m = 0; m < M; ++m)
    {
        const T *const X_row_m = X + m * K;
        T *const Y_row_m = Y + m * N;

        // Initialize Y with bias
        // faster than looping
        std::copy(B, B + N, Y_row_m);

#ifdef INSTRUMENTATION_RUN
        flops += N;
#endif
        for (int k = 0; k < K; ++k)
        {
            const T x_mk_val = X_row_m[k];

            const int W_k_row_start = row_offsets_data[k];
            const int W_k_row_end = row_offsets_data[k + 1];

            for (int nz_idx = W_k_row_start; nz_idx < W_k_row_end; ++nz_idx)
            {
                const int encoded_col_W = encoded_cols_data[nz_idx];
                const int final_col_idx = (encoded_col_W >= 0) ? encoded_col_W : ~encoded_col_W;
                const T val_to_add = (encoded_col_W >= 0) ? x_mk_val : -x_mk_val;

                Y_row_m[final_col_idx] += val_to_add;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }
        }
    }
}

template <typename T>
void ICSC_base(const T *X, const ICSC &W_icsc, const T *B, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_icsc.getDataStructureSize();
#endif

    const int *const col_offsets_data = W_icsc.col_offsets.data();
    const int *const encoded_rows_data = W_icsc.encoded_rows.data();

    for (int m = 0; m < M; ++m)
    {
        const T *const X_row_m = X + m * K;
        T *const Y_row_m = Y + m * N;

        std::copy(B, B + N, Y_row_m);

#ifdef INSTRUMENTATION_RUN
        flops += N;
#endif

        for (int n = 0; n < N; ++n)
        {
            const int W_col_n_start_idx = col_offsets_data[n];
            const int W_col_n_end_idx = col_offsets_data[n + 1];

            for (int nz_idx = W_col_n_start_idx; nz_idx < W_col_n_end_idx; ++nz_idx)
            {
                const int encoded_row_W = encoded_rows_data[nz_idx];
                const int k = (encoded_row_W >= 0) ? encoded_row_W : ~encoded_row_W;
                const T x_mk_val = X_row_m[k];
                const T val_to_add = (encoded_row_W >= 0) ? x_mk_val : -x_mk_val;

                Y_row_m[n] += val_to_add;

#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }
        }
    }
}

template <typename T>
void TCSR_inter(const T *X, const BaseTCSR &W_tcsr, const T *B, T *Y,
                int M, int N, int K)
{

#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_tcsr.getDataStructureSize();
#endif

    const int *const row_start_pos_data = W_tcsr.row_start_pos.data();
    const int *const col_index_pos_data = W_tcsr.col_index_pos.data();
    const int *const row_start_neg_data = W_tcsr.row_start_neg.data();
    const int *const col_index_neg_data = W_tcsr.col_index_neg.data();
    const int W_num_rows = W_tcsr.row_start_pos.size() - 1;

    for (int m = 0; m < M; ++m)
    {
        const T *const X_row_m = X + m * K; // Current row of X
        T *const Y_row_m = Y + m * N;       // Current row of Y

        // Initialize Y with bias
        // faster than looping
        std::copy(B, B + N, Y_row_m);

#ifdef INSTRUMENTATION_RUN
        flops += N;
#endif

        // Iterate over rows of W (r corresponds to columns of X_row_m)
        for (int r = 0; r < W_num_rows; ++r)
        {
            T x_mr_val = X_row_m[r];

            if (x_mr_val == static_cast<T>(0))
            {
                continue;
            }

            // Get pointers/indices for +1 column indices in row r of W
            int p_pos_col_idx = row_start_pos_data[r];
            const int end_pos_col_idx = row_start_pos_data[r + 1];

            // Get pointers/indices for -1 column indices in row r of W
            int p_neg_col_idx = row_start_neg_data[r];
            const int end_neg_col_idx = row_start_neg_data[r + 1];

            int len_pos_cols = end_pos_col_idx - p_pos_col_idx;
            int len_neg_cols = end_neg_col_idx - p_neg_col_idx;

            int common_len_cols = std::min(len_pos_cols, len_neg_cols);

            for (int i = 0; i < common_len_cols - 1; i += 2)
            {
                int c_pos = col_index_pos_data[p_pos_col_idx++];
                Y_row_m[c_pos] += x_mr_val;
                c_pos = col_index_pos_data[p_pos_col_idx++];
                Y_row_m[c_pos] += x_mr_val;

                int c_neg = col_index_neg_data[p_neg_col_idx++];
                Y_row_m[c_neg] -= x_mr_val;
                c_neg = col_index_neg_data[p_neg_col_idx++];
                Y_row_m[c_neg] -= x_mr_val;

#ifdef INSTRUMENTATION_RUN
                flops += 4;
#endif
            }

            while (p_pos_col_idx < end_pos_col_idx)
            {
                const int c_pos = col_index_pos_data[p_pos_col_idx++];
                Y_row_m[c_pos] += x_mr_val;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            while (p_neg_col_idx < end_neg_col_idx)
            {
                const int c_neg = col_index_neg_data[p_neg_col_idx++];
                Y_row_m[c_neg] -= x_mr_val;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }
        }
    }
}

template <typename T>
void TCSC_inter(T *X, const BaseTCSC &W_csc, T *b, T *Y, int M, int N, int K)
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
        const T *X_row_m = X + m * K;
        for (int n = 0; n < N; n++)
        {
            T y_val = 0;

            int p_pos_k_idx = col_start_pos[n];
            const int end_pos_k_idx = col_start_pos[n + 1];

            int p_neg_k_idx = col_start_neg[n];
            const int end_neg_k_idx = col_start_neg[n + 1];

            int len_pos_entries = end_pos_k_idx - p_pos_k_idx;
            int len_neg_entries = end_neg_k_idx - p_neg_k_idx;

            int common_len_entries = std::min(len_pos_entries, len_neg_entries);

            for (int i = 0; i + 1 < common_len_entries; i += 2)
            {
                T x_val_p1 = X_row_m[row_index_pos[p_pos_k_idx++]];
                y_val += x_val_p1;
                T x_val_p2 = X_row_m[row_index_pos[p_pos_k_idx++]];
                y_val += x_val_p2;

                T x_val_n1 = X_row_m[row_index_neg[p_neg_k_idx++]];
                y_val -= x_val_n1;
                T x_val_n2 = X_row_m[row_index_neg[p_neg_k_idx++]];
                y_val -= x_val_n2;

#ifdef INSTRUMENTATION_RUN
                flops += 4; // 2 additions, 2 subtractions
#endif
            }

            while (p_pos_k_idx < end_pos_k_idx)
            {
                T x_val = X_row_m[row_index_pos[p_pos_k_idx++]];
                y_val += x_val;

#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            while (p_neg_k_idx < end_neg_k_idx)
            {
                T x_val = X_row_m[row_index_neg[p_neg_k_idx++]];
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