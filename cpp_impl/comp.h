#ifndef COMP_H
#define COMP_H

#include "common.h"
#include <iostream>
#include <arm_neon.h>

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

// TCSC

template <typename T>
void BaseTCSC(T *X, const TCSC &W_csc, T *b, T *Y, int M, int N, int K)
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

template <typename T, int UNROLL_FACTOR>
void UnrolledModifiedTCSC(T *X, const TCSC &W_csc, T *b, T *Y, int M, int N, int K)
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
            T y_acc[UNROLL_FACTOR];
            for (int u = 0; u < UNROLL_FACTOR; ++u)
            {
                y_acc[u] = T(0);
            }

            int p_pos_k_idx = col_start_pos[n];
            const int end_pos_k_idx = col_start_pos[n + 1];
            int p_neg_k_idx = col_start_neg[n];
            const int end_neg_k_idx = col_start_neg[n + 1];

            while (p_pos_k_idx + UNROLL_FACTOR <= end_pos_k_idx &&
                   p_neg_k_idx + UNROLL_FACTOR <= end_neg_k_idx)
            {
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {

                    y_acc[u] += X_row_m[row_index_pos[p_pos_k_idx + u]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
                p_pos_k_idx += UNROLL_FACTOR;

                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
                    y_acc[u] -= X_row_m[row_index_neg[p_neg_k_idx + u]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
                p_neg_k_idx += UNROLL_FACTOR;
            }

            while (p_pos_k_idx + UNROLL_FACTOR <= end_pos_k_idx)
            {
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
                    y_acc[u] += X_row_m[row_index_pos[p_pos_k_idx + u]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
                p_pos_k_idx += UNROLL_FACTOR;
            }

            while (p_pos_k_idx < end_pos_k_idx)
            {
                y_acc[0] += X_row_m[row_index_pos[p_pos_k_idx++]];
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            while (p_neg_k_idx + UNROLL_FACTOR <= end_neg_k_idx)
            {
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
                    y_acc[u] -= X_row_m[row_index_neg[p_neg_k_idx + u]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
                p_neg_k_idx += UNROLL_FACTOR;
            }

            while (p_neg_k_idx < end_neg_k_idx)
            {
                y_acc[0] -= X_row_m[row_index_neg[p_neg_k_idx++]];
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            T y_final = T(0);
            for (int u = 0; u < UNROLL_FACTOR; u++)
            {
                y_final += y_acc[u];
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            Y[m * N + n] = y_final + b[n];
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif
        }
    }
}

template <typename T, int UNROLL_FACTOR>
void UnrolledTCSC(T *X, const TCSC &W_csc, T *b, T *Y, int M, int N, int K)
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
void BaseInterleavedTCSC(T *X, const InterleavedTCSC &W_csc, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif

    const int *indices_data = W_csc.all_indices.data();
    const int *segment_ptr_data = W_csc.col_segment_ptr.data();

    for (int m = 0; m < M; ++m)
    {
        const T *X_row_m = X + m * K;
        for (int n = 0; n < N; ++n)
        {
            T y_val = 0;

            int pn_start_idx = segment_ptr_data[3 * n + 0];
            int rem_pos_start_idx = segment_ptr_data[3 * n + 1];
            int rem_neg_start_idx = segment_ptr_data[3 * n + 2];
            int next_col_start_idx = segment_ptr_data[3 * n + 3];

            // change +4 to +8 for groups of 4 or + 2 for groups of 1
            for (int k_ptr = pn_start_idx; k_ptr < rem_pos_start_idx; k_ptr += 8)
            {
                y_val += X_row_m[indices_data[k_ptr]];
                y_val += X_row_m[indices_data[k_ptr + 1]];
                y_val += X_row_m[indices_data[k_ptr + 2]];
                y_val += X_row_m[indices_data[k_ptr + 3]];

                // copy and increment index
                y_val -= X_row_m[indices_data[k_ptr + 4]];
                y_val -= X_row_m[indices_data[k_ptr + 5]];
                y_val -= X_row_m[indices_data[k_ptr + 6]];
                y_val -= X_row_m[indices_data[k_ptr + 7]];

                // make sure to change the flops
#ifdef INSTRUMENTATION_RUN
                flops += 8;
#endif
            }

            for (int k_ptr = rem_pos_start_idx; k_ptr < rem_neg_start_idx; ++k_ptr)
            {
                y_val += X_row_m[indices_data[k_ptr]];
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            for (int k_ptr = rem_neg_start_idx; k_ptr < next_col_start_idx; ++k_ptr)
            {
                y_val -= X_row_m[indices_data[k_ptr]];
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
void UnrolledInterleavedTCSC(T *X, const InterleavedTCSC &W_csc, T *b, T *Y, int M, int N, int K)

{

#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif

    const int *indices_data = W_csc.all_indices.data();
    const int *segment_ptr_data = W_csc.col_segment_ptr.data();

    const int BLOCK_SIZE = 8;

    for (int m = 0; m < M; ++m)
    {
        const T *X_row_m = X + m * K;
        for (int n = 0; n < N; ++n)
        {
            T y_val = 0;

            int pn_start_idx = segment_ptr_data[3 * n + 0];
            int rem_pos_start_idx = segment_ptr_data[3 * n + 1];
            int rem_neg_start_idx = segment_ptr_data[3 * n + 2];
            int next_col_start_idx = segment_ptr_data[3 * n + 3];

            int k_ptr;

            // Main unrolled loop
            if (pn_start_idx < rem_pos_start_idx)
            {
                T y_val_accumulators[UNROLL_FACTOR];
                for (int i = 0; i < UNROLL_FACTOR; ++i)
                {
                    y_val_accumulators[i] = 0;
                }

                int num_elements_in_segment1 = rem_pos_start_idx - pn_start_idx;
                int num_blocks = num_elements_in_segment1 / BLOCK_SIZE;
                int num_unrolled_iterations = num_blocks / UNROLL_FACTOR;

                int unrolled_part_limit_loop1 = pn_start_idx + num_unrolled_iterations * UNROLL_FACTOR * BLOCK_SIZE;

                for (k_ptr = pn_start_idx; k_ptr < unrolled_part_limit_loop1; k_ptr += UNROLL_FACTOR * BLOCK_SIZE)
                {
                    for (int u = 0; u < UNROLL_FACTOR; ++u)
                    {
                        int current_block_start_idx = k_ptr + u * BLOCK_SIZE;

                        T val0 = X_row_m[indices_data[current_block_start_idx + 0]];
                        T val1 = X_row_m[indices_data[current_block_start_idx + 1]];
                        T val2 = X_row_m[indices_data[current_block_start_idx + 2]];
                        T val3 = X_row_m[indices_data[current_block_start_idx + 3]];
                        T val4 = X_row_m[indices_data[current_block_start_idx + 4]];
                        T val5 = X_row_m[indices_data[current_block_start_idx + 5]];
                        T val6 = X_row_m[indices_data[current_block_start_idx + 6]];
                        T val7 = X_row_m[indices_data[current_block_start_idx + 7]];

                        T p_sum1 = val0 + val1;
                        T p_sum2 = val2 + val3;
                        T n_sum1 = val4 + val5;
                        T n_sum2 = val6 + val7;

                        T total_pos = p_sum1 + p_sum2;
                        T total_neg = n_sum1 + n_sum2;

                        y_val_accumulators[u] += total_pos - total_neg;

#ifdef INSTRUMENTATION_RUN
                        flops += BLOCK_SIZE;
#endif
                    }
                }

                for (int u = 0; u < UNROLL_FACTOR; ++u)
                {
                    y_val += y_val_accumulators[u];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                for (; k_ptr < rem_pos_start_idx; k_ptr += BLOCK_SIZE)
                {
                    // for cleanup
                    T val0 = X_row_m[indices_data[k_ptr + 0]];
                    T val1 = X_row_m[indices_data[k_ptr + 1]];
                    T val2 = X_row_m[indices_data[k_ptr + 2]];
                    T val3 = X_row_m[indices_data[k_ptr + 3]];
                    T val4 = X_row_m[indices_data[k_ptr + 4]];
                    T val5 = X_row_m[indices_data[k_ptr + 5]];
                    T val6 = X_row_m[indices_data[k_ptr + 6]];
                    T val7 = X_row_m[indices_data[k_ptr + 7]];

                    T p_sum1 = val0 + val1;
                    T p_sum2 = val2 + val3;
                    T n_sum1 = val4 + val5;
                    T n_sum2 = val6 + val7;

                    T total_pos = p_sum1 + p_sum2;
                    T total_neg = n_sum1 + n_sum2;

                    y_val += total_pos - total_neg;
#ifdef INSTRUMENTATION_RUN
                    flops += BLOCK_SIZE;
#endif
                }
            }

            // Remaining positive
            if (rem_pos_start_idx < rem_neg_start_idx)
            {
                for (k_ptr = rem_pos_start_idx; k_ptr < rem_neg_start_idx; ++k_ptr)
                {
                    y_val += X_row_m[indices_data[k_ptr]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
            }

            // Remaining negative
            if (rem_neg_start_idx < next_col_start_idx)
            {
                for (k_ptr = rem_neg_start_idx; k_ptr < next_col_start_idx; ++k_ptr)
                {
                    y_val -= X_row_m[indices_data[k_ptr]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
            }

            Y[m * N + n] = y_val + b[n];
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif
        }
    }
}

// TCSR

template <typename T>
void BaseTCSR(T *X, const TCSR &W_csr, T *b, T *Y, int M, int N, int K)
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
void UnrolledBaseTCSR(T *X, const TCSR &W_csr, T *b, T *Y, int M, int N, int K)
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

// BlockedTCSC

template <typename T, int B>
void BaseBlockedTCSC(T *X, const BlockedTCSC<B> &W_csc, T *b, T *Y, int M, int N, int K)
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
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
                for (int k = col_start_neg[n]; k < col_start_neg[n + 1]; k++)
                {
                    y -= X[m * K + row_index_neg[k]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
                Y[m * N + n % N] += y;
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

template <typename T, int B, int UNROLL_FACTOR>
void UnrolledBlockedTCSC(T *X, const BlockedTCSC<B> &W_csc, T *b, T *Y, int M, int N, int K)
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
        const T *X_row_m = X + m * K;

        // Initialize Y row with bias (only once)
        for (int n = 0; n < N; n++)
        {
            Y[m * N + n] = b[n];
        }

        // Temporary accumulator for each output column to avoid WAW hazards
        T y_temp[N];
        for (int n = 0; n < N; n++)
        {
            y_temp[n] = T(0);
        }

        // Process each block of K
        for (int k_block = 0; k_block < K / B; k_block++)
        {
            for (int n = k_block * N; n < k_block * N + N; n++)
            {
                int output_col = n % N; // Which output column we're contributing to

                // Multiple accumulators to break dependency chains
                T y_acc[UNROLL_FACTOR];
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
                    y_acc[u] = T(0);
                }

                // Process positive values with unrolling
                int k_pos = col_start_pos[n];
                const int end_pos = col_start_pos[n + 1];

                // Main unrolled loop for positive values
                for (; k_pos + UNROLL_FACTOR <= end_pos; k_pos += UNROLL_FACTOR)
                {
                    for (int u = 0; u < UNROLL_FACTOR; u++)
                    {
                        y_acc[u] += X_row_m[row_index_pos[k_pos + u]];
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                    }
                }

                // Cleanup remaining positive values
                for (; k_pos < end_pos; k_pos++)
                {
                    y_acc[0] += X_row_m[row_index_pos[k_pos]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                // Process negative values with unrolling
                int k_neg = col_start_neg[n];
                const int end_neg = col_start_neg[n + 1];

                // Main unrolled loop for negative values
                for (; k_neg + UNROLL_FACTOR <= end_neg; k_neg += UNROLL_FACTOR)
                {
                    for (int u = 0; u < UNROLL_FACTOR; u++)
                    {
                        y_acc[u] -= X_row_m[row_index_neg[k_neg + u]];
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                    }
                }

                // Cleanup remaining negative values
                for (; k_neg < end_neg; k_neg++)
                {
                    y_acc[0] -= X_row_m[row_index_neg[k_neg]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                // Combine all accumulators and add to temporary
                T y_final = T(0);
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
                    y_final += y_acc[u];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                y_temp[output_col] += y_final; // Accumulate into temp (no WAW hazard)
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }
        }

        // Final write to Y (eliminates WAW hazards)
        for (int n = 0; n < N; n++)
        {
            Y[m * N + n] += y_temp[n];
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif
        }
    }
}

template <typename T, int B>
void BaseInterleavedBlockedTCSC(T *X, const InterleavedBlockedTCSC<B> &W_csc, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif

    const int *indices_data = W_csc.all_indices.data();
    const int *segment_ptr_data = W_csc.col_segment_ptr.data();

    int num_blocks = K / B;

    // Process each row of X
    for (int m = 0; m < M; m++)
    {
        for (int j = 0; j < N; j++)
            Y[m * N + j] = b[j];
        // Process each column-block of X
        T *X_row_m = X + m * K;
        for (int k_block = 0; k_block < num_blocks; k_block++)
        {
            for (int j = 0; j < N; j++)
            {
                int index = k_block * N + j;
                T sumPos = 0;
                T sumNeg = 0;

                int pn_start_idx = segment_ptr_data[3 * index + 0];
                int rem_pos_start_idx = segment_ptr_data[3 * index + 1];
                int rem_neg_start_idx = segment_ptr_data[3 * index + 2];
                int next_col_start_idx = segment_ptr_data[3 * index + 3];

                // change +4 to +8 for groups of 4 or + 2 for groups of 1
                for (int k_ptr = pn_start_idx; k_ptr < rem_pos_start_idx; k_ptr += 2)
                {
                    sumPos += X_row_m[indices_data[k_ptr]];
                    sumNeg -= X_row_m[indices_data[k_ptr + 1]];
#ifdef INSTRUMENTATION_RUN
                    flops += 2;
#endif
                }

                int k_ptr = rem_pos_start_idx;

                while (k_ptr < rem_neg_start_idx)
                {
                    sumPos += X_row_m[indices_data[k_ptr++]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                while (k_ptr < next_col_start_idx)
                {
                    sumNeg -= X_row_m[indices_data[k_ptr++]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                T y_val = sumPos + sumNeg;
                Y[m * N + j] += y_val;
            }
        }
    }
}

// UNROLL_FACTOR must divide N
template <typename T, int B, int UNROLL_FACTOR>
void UnrolledInterleavedBlockedTCSC(T *X, const InterleavedBlockedTCSC<B> &W_csc, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    constexpr int UNROLL_FACTOR_HALF = UNROLL_FACTOR / 2;
    const int *indices_data = W_csc.all_indices.data();
    const int *segment_ptr_data = W_csc.col_segment_ptr.data();

    int num_blocks = K / B;

    // Process each row of X
    for (int m = 0; m < M; m++)
    {
        for (int j = 0; j < N; j++)
            Y[m * N + j] = b[j];

        // Process each column-block of X
        T *X_row_m = X + m * K;
        for (int k_block = 0; k_block < num_blocks; k_block++)
        {
            for (int j = 0; j < N; j++)
            {
                T y_acc[UNROLL_FACTOR];
                for (int u = 0; u < UNROLL_FACTOR; ++u)
                    y_acc[u] = 0;

                int base_index = k_block * N + j;

                int pn_start_idx = segment_ptr_data[3 * base_index + 0];
                int rem_pos_start_idx = segment_ptr_data[3 * base_index + 1];
                int rem_neg_start_idx = segment_ptr_data[3 * base_index + 2];
                int next_col_start_idx = segment_ptr_data[3 * base_index + 3];

                for (int k_ptr = pn_start_idx; k_ptr < rem_pos_start_idx; k_ptr += UNROLL_FACTOR)
                {
                    for (int i = 0, p = UNROLL_FACTOR_HALF; i < UNROLL_FACTOR_HALF; i++, p++)
                    {
                        y_acc[i] += X_row_m[indices_data[k_ptr + i]];
                        y_acc[p] -= X_row_m[indices_data[k_ptr + p]];
#ifdef INSTRUMENTATION_RUN
                        flops += 2;
#endif
                    }
                }
                int k_ptr = rem_pos_start_idx;

                while (k_ptr < rem_neg_start_idx)
                {
                    y_acc[0] += X_row_m[indices_data[k_ptr++]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                while (k_ptr < next_col_start_idx)
                {
                    y_acc[1] -= X_row_m[indices_data[k_ptr++]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }

                // Here I tried to unroll the previously serialized addition
                constexpr int mid = UNROLL_FACTOR_HALF;
                T y_acc0 = 0, y_acc1 = 0;
                for (int u = 0, p = mid; u < UNROLL_FACTOR_HALF; u++, p++)
                {
                    y_acc0 += y_acc[u];
                    y_acc1 += y_acc[p];
#ifdef INSTRUMENTATION_RUN
                    flops += 2;
#endif
                }

                Y[m * N + j] += y_acc0 + y_acc1;
            }
        }
    }
}

template <typename T>
void NeonInterleavedTCSC(T *X, const InterleavedTCSC &W_csc, T *b, T *Y, int M, int N, int K)
{
    const int *indices_data = W_csc.all_indices.data();
    const int *segment_ptr_data = W_csc.col_segment_ptr.data();

    // Create a sign vector for addition and subtraction
    const float32x4_t signs = {1.0f, 1.0f, -1.0f, -1.0f};

    for (int m = 0; m < M; ++m)
    {
        const T *X_row_m = X + m * K;
        for (int n = 0; n < N; n += 4)
        {
            T y_val0 = 0;
            T y_val1 = 0;
            T y_val2 = 0;
            T y_val3 = 0;
            { // n = 0
                // Initialize an accumulator vector to zeros
                float32x4_t y_val_vec = vdupq_n_f32(0.0f);

                // Load counters from ds
                int32x4_t indices_vec = vld1q_s32(segment_ptr_data + 3 * n);
                int pn_start_idx = vgetq_lane_s32(indices_vec, 0);
                int rem_pos_start_idx = vgetq_lane_s32(indices_vec, 1);
                int rem_neg_start_idx = vgetq_lane_s32(indices_vec, 2);
                int next_col_start_idx = vgetq_lane_s32(indices_vec, 3);

                int k_ptr = pn_start_idx;
                for (; k_ptr + 3 < rem_pos_start_idx; k_ptr += 4)
                {
                    // Load indices
                    int32x4_t indices_vec = vld1q_s32(indices_data + k_ptr);

                    // Gather values from X_row_m using the indices
                    float32x4_t x_vals = {
                        X_row_m[vgetq_lane_s32(indices_vec, 0)],
                        X_row_m[vgetq_lane_s32(indices_vec, 1)],
                        X_row_m[vgetq_lane_s32(indices_vec, 2)],
                        X_row_m[vgetq_lane_s32(indices_vec, 3)]};

                    // Perform fused multiply-add: y_val_vec += x_vals * signs
                    y_val_vec = vmlaq_f32(y_val_vec, x_vals, signs);
                }

                // Horizontally add the elements of the accumulator vector
                y_val0 += vaddvq_f32(y_val_vec);

                // Handle the remaining elements (positive contributions)
                for (k_ptr = rem_pos_start_idx; k_ptr < rem_neg_start_idx; ++k_ptr)
                {
                    y_val0 += X_row_m[indices_data[k_ptr]];
                }

                // Handle the remaining elements (negative contributions)
                for (k_ptr = rem_neg_start_idx; k_ptr < next_col_start_idx; ++k_ptr)
                {
                    y_val0 -= X_row_m[indices_data[k_ptr]];
                }
            }

            { // n = 1
                // Initialize an accumulator vector to zeros
                float32x4_t y_val_vec = vdupq_n_f32(0.0f);

                // Load counters from ds
                int32x4_t indices_vec = vld1q_s32(segment_ptr_data + 3 * (n + 1));
                int pn_start_idx = vgetq_lane_s32(indices_vec, 0);
                int rem_pos_start_idx = vgetq_lane_s32(indices_vec, 1);
                int rem_neg_start_idx = vgetq_lane_s32(indices_vec, 2);
                int next_col_start_idx = vgetq_lane_s32(indices_vec, 3);

                int k_ptr = pn_start_idx;
                for (; k_ptr + 3 < rem_pos_start_idx; k_ptr += 4)
                {
                    // Load indices
                    int32x4_t indices_vec = vld1q_s32(indices_data + k_ptr);

                    // Gather values from X_row_m using the indices
                    float32x4_t x_vals = {
                        X_row_m[vgetq_lane_s32(indices_vec, 0)],
                        X_row_m[vgetq_lane_s32(indices_vec, 1)],
                        X_row_m[vgetq_lane_s32(indices_vec, 2)],
                        X_row_m[vgetq_lane_s32(indices_vec, 3)]};

                    // Perform fused multiply-add: y_val_vec += x_vals * signs
                    y_val_vec = vmlaq_f32(y_val_vec, x_vals, signs);
                }

                // Horizontally add the elements of the accumulator vector
                y_val1 += vaddvq_f32(y_val_vec);

                // Handle the remaining elements (positive contributions)
                for (k_ptr = rem_pos_start_idx; k_ptr < rem_neg_start_idx; ++k_ptr)
                {
                    y_val1 += X_row_m[indices_data[k_ptr]];
                }

                // Handle the remaining elements (negative contributions)
                for (k_ptr = rem_neg_start_idx; k_ptr < next_col_start_idx; ++k_ptr)
                {
                    y_val1 -= X_row_m[indices_data[k_ptr]];
                }
            }

            { // n = 2
                // Initialize an accumulator vector to zeros
                float32x4_t y_val_vec = vdupq_n_f32(0.0f);

                // Load counters from ds
                int32x4_t indices_vec = vld1q_s32(segment_ptr_data + 3 * (n + 2));
                int pn_start_idx = vgetq_lane_s32(indices_vec, 0);
                int rem_pos_start_idx = vgetq_lane_s32(indices_vec, 1);
                int rem_neg_start_idx = vgetq_lane_s32(indices_vec, 2);
                int next_col_start_idx = vgetq_lane_s32(indices_vec, 3);

                int k_ptr = pn_start_idx;
                for (; k_ptr + 3 < rem_pos_start_idx; k_ptr += 4)
                {
                    // Load indices
                    int32x4_t indices_vec = vld1q_s32(indices_data + k_ptr);

                    // Gather values from X_row_m using the indices
                    float32x4_t x_vals = {
                        X_row_m[vgetq_lane_s32(indices_vec, 0)],
                        X_row_m[vgetq_lane_s32(indices_vec, 1)],
                        X_row_m[vgetq_lane_s32(indices_vec, 2)],
                        X_row_m[vgetq_lane_s32(indices_vec, 3)]};

                    // Perform fused multiply-add: y_val_vec += x_vals * signs
                    y_val_vec = vmlaq_f32(y_val_vec, x_vals, signs);
                }

                // Horizontally add the elements of the accumulator vector
                y_val2 += vaddvq_f32(y_val_vec);

                // Handle the remaining elements (positive contributions)
                for (k_ptr = rem_pos_start_idx; k_ptr < rem_neg_start_idx; ++k_ptr)
                {
                    y_val2 += X_row_m[indices_data[k_ptr]];
                }

                // Handle the remaining elements (negative contributions)
                for (k_ptr = rem_neg_start_idx; k_ptr < next_col_start_idx; ++k_ptr)
                {
                    y_val2 -= X_row_m[indices_data[k_ptr]];
                }
            }
            { // n = 3
                // Initialize an accumulator vector to zeros
                float32x4_t y_val_vec = vdupq_n_f32(0.0f);

                // Load counters from ds
                int32x4_t indices_vec = vld1q_s32(segment_ptr_data + 3 * (n + 3));
                int pn_start_idx = vgetq_lane_s32(indices_vec, 0);
                int rem_pos_start_idx = vgetq_lane_s32(indices_vec, 1);
                int rem_neg_start_idx = vgetq_lane_s32(indices_vec, 2);
                int next_col_start_idx = vgetq_lane_s32(indices_vec, 3);

                int k_ptr = pn_start_idx;
                for (; k_ptr + 3 < rem_pos_start_idx; k_ptr += 4)
                {
                    // Load indices
                    int32x4_t indices_vec = vld1q_s32(indices_data + k_ptr);

                    // Gather values from X_row_m using the indices
                    float32x4_t x_vals = {
                        X_row_m[vgetq_lane_s32(indices_vec, 0)],
                        X_row_m[vgetq_lane_s32(indices_vec, 1)],
                        X_row_m[vgetq_lane_s32(indices_vec, 2)],
                        X_row_m[vgetq_lane_s32(indices_vec, 3)]};

                    // Perform fused multiply-add: y_val_vec += x_vals * signs
                    y_val_vec = vmlaq_f32(y_val_vec, x_vals, signs);
                }

                // Horizontally add the elements of the accumulator vector
                y_val3 += vaddvq_f32(y_val_vec);

                // Handle the remaining elements (positive contributions)
                for (k_ptr = rem_pos_start_idx; k_ptr < rem_neg_start_idx; ++k_ptr)
                {
                    y_val3 += X_row_m[indices_data[k_ptr]];
                }

                // Handle the remaining elements (negative contributions)
                for (k_ptr = rem_neg_start_idx; k_ptr < next_col_start_idx; ++k_ptr)
                {
                    y_val3 -= X_row_m[indices_data[k_ptr]];
                }
            }

            float32x4_t y_res = {y_val0, y_val1, y_val2, y_val3};
            float32x4_t b_vals = vld1q_f32(b + n);
            float32x4_t store_in_y = vaddq_f32(y_res, b_vals);
            vst1q_f32(Y + m * N + n, store_in_y);
        }
    }
}

#endif