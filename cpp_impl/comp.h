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
void UnrolledSimultaneousTCSC(T *X, const TCSC &W_csc, T *b, T *Y, int M, int N, int K)
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

template <typename T, int UNROLL_FACTOR>
void UnrolledTCSR(T *X, const TCSR &W_csr, T *b, T *Y, int M, int N, int K)
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
void NeonTCSCHorizontalSimple(float *X, const TCSC &W_csc, float *b, float *Y, int M, int N, int K)
{
    const int *col_start_pos = W_csc.col_start_pos.data();
    const int *col_start_neg = W_csc.col_start_neg.data();
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();

    for (int m = 0; m < M; m++)
    {
        const float *X_row_m = X + m * K;

        for (int n = 0; n < N; n += 4)
        {
            float32x4_t y_vec0 = vdupq_n_f32(0.0f);
            float32x4_t y_vec1 = vdupq_n_f32(0.0f);
            float32x4_t y_vec2 = vdupq_n_f32(0.0f);
            float32x4_t y_vec3 = vdupq_n_f32(0.0f);

            // Process positive values
            int k0 = col_start_pos[n];
            int k1 = col_start_pos[n + 1];
            int k2 = col_start_pos[n + 2];
            int k3 = col_start_pos[n + 3];

            int k0_end = col_start_pos[n + 1];
            for (; k0 + 3 < k0_end; k0 += 4)
            {
                float32x4_t x_vals = {X_row_m[row_index_pos[k0]], X_row_m[row_index_pos[k0 + 1]], X_row_m[row_index_pos[k0 + 2]], X_row_m[row_index_pos[k0 + 3]]};
                y_vec0 = vaddq_f32(y_vec0, x_vals);
            }
            int k1_end = col_start_pos[n + 2];
            for (; k1 + 3 < k1_end; k1 += 4)
            {
                float32x4_t x_vals = {X_row_m[row_index_pos[k1]], X_row_m[row_index_pos[k1 + 1]], X_row_m[row_index_pos[k1 + 2]], X_row_m[row_index_pos[k1 + 3]]};
                y_vec1 = vaddq_f32(y_vec1, x_vals);
            }
            int k2_end = col_start_pos[n + 3];
            for (; k2 + 3 < k2_end; k2 += 4)
            {
                float32x4_t x_vals = {X_row_m[row_index_pos[k2]], X_row_m[row_index_pos[k2 + 1]], X_row_m[row_index_pos[k2 + 2]], X_row_m[row_index_pos[k2 + 3]]};
                y_vec2 = vaddq_f32(y_vec2, x_vals);
            }
            int k3_end = col_start_pos[n + 4];
            for (; k3 + 3 < k3_end; k3 += 4)
            {
                float32x4_t x_vals = {X_row_m[row_index_pos[k3]], X_row_m[row_index_pos[k3 + 1]], X_row_m[row_index_pos[k3 + 2]], X_row_m[row_index_pos[k3 + 3]]};
                y_vec3 = vaddq_f32(y_vec3, x_vals);
            }

            // Process negative values
            int l0 = col_start_neg[n];
            int l1 = col_start_neg[n + 1];
            int l2 = col_start_neg[n + 2];
            int l3 = col_start_neg[n + 3];

            int l0_end = col_start_neg[n + 1];
            for (; l0 + 3 < l0_end; l0 += 4)
            {
                float32x4_t x_vals = {X_row_m[row_index_neg[l0]], X_row_m[row_index_neg[l0 + 1]], X_row_m[row_index_neg[l0 + 2]], X_row_m[row_index_neg[l0 + 3]]};
                y_vec0 = vsubq_f32(y_vec0, x_vals);
            }
            int l1_end = col_start_neg[n + 2];
            for (; l1 + 3 < l1_end; l1 += 4)
            {
                float32x4_t x_vals = {X_row_m[row_index_neg[l1]], X_row_m[row_index_neg[l1 + 1]], X_row_m[row_index_neg[l1 + 2]], X_row_m[row_index_neg[l1 + 3]]};
                y_vec1 = vsubq_f32(y_vec1, x_vals);
            }
            int l2_end = col_start_neg[n + 3];
            for (; l2 + 3 < l2_end; l2 += 4)
            {
                float32x4_t x_vals = {X_row_m[row_index_neg[l2]], X_row_m[row_index_neg[l2 + 1]], X_row_m[row_index_neg[l2 + 2]], X_row_m[row_index_neg[l2 + 3]]};
                y_vec2 = vsubq_f32(y_vec2, x_vals);
            }
            int l3_end = col_start_neg[n + 4];
            for (; l3 + 3 < l3_end; l3 += 4)
            {
                float32x4_t x_vals = {X_row_m[row_index_neg[l3]], X_row_m[row_index_neg[l3 + 1]], X_row_m[row_index_neg[l3 + 2]], X_row_m[row_index_neg[l3 + 3]]};
                y_vec3 = vsubq_f32(y_vec3, x_vals);
            }

            // Horizontal add
            float y0 = vaddvq_f32(y_vec0);
            float y1 = vaddvq_f32(y_vec1);
            float y2 = vaddvq_f32(y_vec2);
            float y3 = vaddvq_f32(y_vec3);

            // Remainder Loops
            for (; k0 < k0_end; k0++)
            {
                y0 += X_row_m[row_index_pos[k0]];
            }
            for (; l0 < l0_end; l0++)
            {
                y0 -= X_row_m[row_index_neg[l0]];
            }
            for (; k1 < k1_end; k1++)
            {
                y1 += X_row_m[row_index_pos[k1]];
            }
            for (; l1 < l1_end; l1++)
            {
                y1 -= X_row_m[row_index_neg[l1]];
            }
            for (; k2 < k2_end; k2++)
            {
                y2 += X_row_m[row_index_pos[k2]];
            }
            for (; l2 < l2_end; l2++)
            {
                y2 -= X_row_m[row_index_neg[l2]];
            }
            for (; k3 < k3_end; k3++)
            {
                y3 += X_row_m[row_index_pos[k3]];
            }
            for (; l3 < l3_end; l3++)
            {
                y3 -= X_row_m[row_index_neg[l3]];
            }

            // Store final results
            Y[m * N + n] = y0 + b[n];
            Y[m * N + n + 1] = y1 + b[n + 1];
            Y[m * N + n + 2] = y2 + b[n + 2];
            Y[m * N + n + 3] = y3 + b[n + 3];
            // float32x4_t y_res = { y0, y1, y2, y3 };
            // float32x4_t b_vals = vld1q_f32(b + n);
            // float32x4_t store_in_y = vaddq_f32(y_res, b_vals);
            // vst1q_f32(Y + m * N + n, store_in_y);
        }
    }
}

template <typename T>
void NeonTCSCVertical(float *X, const VectorTCSC &W_csc, float *b, float *Y, int M, int N, int K)
{
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();
    const int *cap_every_four = W_csc.cap_every_four.data();

    for (int m = 0; m < M; m++)
    {
        float *X_row_m = X + m * K;
        X_row_m[-1] = 0;

        int cap_idx = 0;
        int k3_end = 0;
        for (int n = 0; n < N; n += 4)
        {
            float32x4_t y_vec0 = vdupq_n_f32(0.0f);
            float32x4_t y_vec1 = vdupq_n_f32(0.0f);

            int cap = cap_every_four[cap_idx++];
            int k0 = k3_end;
            int k0_end = k0 + cap;
            int k1 = k0 + cap;
            int k2 = k1 + cap;
            int k3 = k2 + cap;
            k3_end = k3 + cap;

            for (int i = 0; i < cap; i += 4)
            {
                float32x4_t x_vals000 = {X_row_m[row_index_pos[k0]],
                                         X_row_m[row_index_pos[k1]],
                                         X_row_m[row_index_pos[k2]],
                                         X_row_m[row_index_pos[k3]]};
                float32x4_t x_vals001 = {X_row_m[row_index_neg[k0]],
                                         X_row_m[row_index_neg[k1]],
                                         X_row_m[row_index_neg[k2]],
                                         X_row_m[row_index_neg[k3]]};
                float32x4_t x_vals010 = {X_row_m[row_index_pos[k0 + 1]],
                                         X_row_m[row_index_pos[k1 + 1]],
                                         X_row_m[row_index_pos[k2 + 1]],
                                         X_row_m[row_index_pos[k3 + 1]]};
                float32x4_t x_vals011 = {X_row_m[row_index_neg[k0 + 1]],
                                         X_row_m[row_index_neg[k1 + 1]],
                                         X_row_m[row_index_neg[k2 + 1]],
                                         X_row_m[row_index_neg[k3 + 1]]};
                float32x4_t x_vals100 = {X_row_m[row_index_pos[k0 + 2]],
                                         X_row_m[row_index_pos[k1 + 2]],
                                         X_row_m[row_index_pos[k2 + 2]],
                                         X_row_m[row_index_pos[k3 + 2]]};
                float32x4_t x_vals101 = {X_row_m[row_index_neg[k0 + 2]],
                                         X_row_m[row_index_neg[k1 + 2]],
                                         X_row_m[row_index_neg[k2 + 2]],
                                         X_row_m[row_index_neg[k3 + 2]]};
                float32x4_t x_vals110 = {X_row_m[row_index_pos[k0 + 3]],
                                         X_row_m[row_index_pos[k1 + 3]],
                                         X_row_m[row_index_pos[k2 + 3]],
                                         X_row_m[row_index_pos[k3 + 3]]};
                float32x4_t x_vals111 = {X_row_m[row_index_neg[k0 + 3]],
                                         X_row_m[row_index_neg[k1 + 3]],
                                         X_row_m[row_index_neg[k2 + 3]],
                                         X_row_m[row_index_neg[k3 + 3]]};
                k0 += 4;
                k1 += 4;
                k2 += 4;
                k3 += 4;
                float32x4_t vec_tmp0 = vsubq_f32(x_vals000, x_vals001);
                float32x4_t vec_tmp1 = vsubq_f32(x_vals010, x_vals011);
                float32x4_t vec_tmp2 = vsubq_f32(x_vals100, x_vals101);
                float32x4_t vec_tmp3 = vsubq_f32(x_vals110, x_vals111);
                float32x4_t vec_tmp4 = vaddq_f32(vec_tmp0, vec_tmp1);
                float32x4_t vec_tmp5 = vaddq_f32(vec_tmp2, vec_tmp3);
                y_vec0 = vaddq_f32(y_vec0, vec_tmp4);
                y_vec1 = vaddq_f32(y_vec1, vec_tmp5);
            }

            float32x4_t b_vals = vld1q_f32(b + n);
            float32x4_t y_vec_res = vaddq_f32(y_vec0, y_vec1);
            float32x4_t store_in_y = vaddq_f32(y_vec_res, b_vals);
            vst1q_f32(Y + m * N + n, store_in_y);
        }
    }
}

template <typename T>
void NeonTCSCHorizontalAdvanced(float *X, const VectorTCSC &W_csc, float *b, float *Y, int M, int N, int K)
{
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();
    const int *cap_every_four = W_csc.cap_every_four.data();

    for (int m = 0; m < M; m++)
    {
        float *X_row_m = X + m * K;
        X_row_m[-1] = 0;

        int cap_idx = 0;
        int k3_end = 0;
        for (int n = 0; n < N; n += 4)
        {
            float32x4_t y_vec0 = vdupq_n_f32(0.0f);
            float32x4_t y_vec1 = vdupq_n_f32(0.0f);
            float32x4_t y_vec2 = vdupq_n_f32(0.0f);
            float32x4_t y_vec3 = vdupq_n_f32(0.0f);

            int cap = cap_every_four[cap_idx++];
            int k0 = k3_end;
            int k0_end = k0 + cap;
            int k1 = k0 + cap;
            int k2 = k1 + cap;
            int k3 = k2 + cap;
            k3_end = k3 + cap;

            for (int i = k0; i < k0_end; i += 4)
            {
                float32x4_t x_vals000 = {X_row_m[row_index_pos[k0]], X_row_m[row_index_pos[k0 + 1]], X_row_m[row_index_pos[k0 + 2]], X_row_m[row_index_pos[k0 + 3]]};
                float32x4_t x_vals001 = {X_row_m[row_index_neg[k0]], X_row_m[row_index_neg[k0 + 1]], X_row_m[row_index_neg[k0 + 2]], X_row_m[row_index_neg[k0 + 3]]};
                float32x4_t x_vals010 = {X_row_m[row_index_pos[k1]], X_row_m[row_index_pos[k1 + 1]], X_row_m[row_index_pos[k1 + 2]], X_row_m[row_index_pos[k1 + 3]]};
                float32x4_t x_vals011 = {X_row_m[row_index_neg[k1]], X_row_m[row_index_neg[k1 + 1]], X_row_m[row_index_neg[k1 + 2]], X_row_m[row_index_neg[k1 + 3]]};
                float32x4_t x_vals100 = {X_row_m[row_index_pos[k2]], X_row_m[row_index_pos[k2 + 1]], X_row_m[row_index_pos[k2 + 2]], X_row_m[row_index_pos[k2 + 3]]};
                float32x4_t x_vals101 = {X_row_m[row_index_neg[k2]], X_row_m[row_index_neg[k2 + 1]], X_row_m[row_index_neg[k2 + 2]], X_row_m[row_index_neg[k2 + 3]]};
                float32x4_t x_vals110 = {X_row_m[row_index_pos[k3]], X_row_m[row_index_pos[k3 + 1]], X_row_m[row_index_pos[k3 + 2]], X_row_m[row_index_pos[k3 + 3]]};
                float32x4_t x_vals111 = {X_row_m[row_index_neg[k3]], X_row_m[row_index_neg[k3 + 1]], X_row_m[row_index_neg[k3 + 2]], X_row_m[row_index_neg[k3 + 3]]};

                float32x4_t vec_tmp0 = vsubq_f32(x_vals000, x_vals001);
                float32x4_t vec_tmp1 = vsubq_f32(x_vals010, x_vals011);
                float32x4_t vec_tmp2 = vsubq_f32(x_vals100, x_vals101);
                float32x4_t vec_tmp3 = vsubq_f32(x_vals110, x_vals111);

                k0 += 4;
                k1 += 4;
                k2 += 4;
                k3 += 4;

                y_vec0 = vaddq_f32(y_vec0, vec_tmp0);
                y_vec1 = vaddq_f32(y_vec1, vec_tmp1);
                y_vec2 = vaddq_f32(y_vec2, vec_tmp2);
                y_vec3 = vaddq_f32(y_vec3, vec_tmp3);
            }

            // Horizontal add
            float y0 = vaddvq_f32(y_vec0);
            float y1 = vaddvq_f32(y_vec1);
            float y2 = vaddvq_f32(y_vec2);
            float y3 = vaddvq_f32(y_vec3);

            // Store final results
            Y[m * N + n] = y0 + b[n];
            Y[m * N + n + 1] = y1 + b[n + 1];
            Y[m * N + n + 2] = y2 + b[n + 2];
            Y[m * N + n + 3] = y3 + b[n + 3];
            // // Code below doesn't add much in performance smh
            // float32x4_t y_res = { y0, y1, y2, y3 };
            // float32x4_t b_vals = vld1q_f32(b + n);
            // float32x4_t store_in_y = vaddq_f32(y_res, b_vals);
            // vst1q_f32(Y + m * N + n, store_in_y);
        }
    }
}

template <typename T, int K_UNROLL_FACTOR, int M_UNROLL_FACTOR>
void DoubleUnrolledTCSC(T *X, const TCSC &W_csc, T *b, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *col_start_pos = W_csc.col_start_pos.data();
    const int *col_start_neg = W_csc.col_start_neg.data();
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();

    // Main loop with M unrolling
    int m;
    for (m = 0; m <= M - M_UNROLL_FACTOR; m += M_UNROLL_FACTOR)
    {
        for (int n = 0; n < N; n++)
        {
            // Separate accumulators for each unrolled M iteration
            T y_pos[M_UNROLL_FACTOR][K_UNROLL_FACTOR];
            T y_neg[M_UNROLL_FACTOR][K_UNROLL_FACTOR];

            // Initialize accumulators
            for (int v = 0; v < M_UNROLL_FACTOR; v++)
            {
                for (int u = 0; u < K_UNROLL_FACTOR; u++)
                {
                    y_pos[v][u] = T(0);
                    y_neg[v][u] = T(0);
                }
            }

            // Process positive values with double unrolling
            int k_pos_loop = col_start_pos[n];
            const int end_pos = col_start_pos[n + 1];

            for (; k_pos_loop + K_UNROLL_FACTOR <= end_pos; k_pos_loop += K_UNROLL_FACTOR)
            {
                for (int u = 0; u < K_UNROLL_FACTOR; u++)
                {
                    int row_idx = row_index_pos[k_pos_loop + u];
                    for (int v = 0; v < M_UNROLL_FACTOR; v++)
                    {
                        y_pos[v][u] += X[(m + v) * K + row_idx];
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                    }
                }
            }

            // Reduce positive accumulators for each M iteration
            T y_pos_final[M_UNROLL_FACTOR];
            for (int v = 0; v < M_UNROLL_FACTOR; v++)
            {
                y_pos_final[v] = T(0);
                for (int u = 0; u < K_UNROLL_FACTOR; u++)
                {
                    y_pos_final[v] += y_pos[v][u];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
            }

            // Handle remaining positive values
            for (; k_pos_loop < end_pos; k_pos_loop++)
            {
                int row_idx = row_index_pos[k_pos_loop];
                for (int v = 0; v < M_UNROLL_FACTOR; v++)
                {
                    y_pos_final[v] += X[(m + v) * K + row_idx];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
            }

            // Process negative values with double unrolling
            int k_neg_loop = col_start_neg[n];
            const int end_neg = col_start_neg[n + 1];

            for (; k_neg_loop + K_UNROLL_FACTOR <= end_neg; k_neg_loop += K_UNROLL_FACTOR)
            {
                for (int u = 0; u < K_UNROLL_FACTOR; u++)
                {
                    int row_idx = row_index_neg[k_neg_loop + u];
                    for (int v = 0; v < M_UNROLL_FACTOR; v++)
                    {
                        y_neg[v][u] += X[(m + v) * K + row_idx];
#ifdef INSTRUMENTATION_RUN
                        flops++;
#endif
                    }
                }
            }

            // Reduce negative accumulators for each M iteration
            T y_neg_final[M_UNROLL_FACTOR];
            for (int v = 0; v < M_UNROLL_FACTOR; v++)
            {
                y_neg_final[v] = T(0);
                for (int u = 0; u < K_UNROLL_FACTOR; u++)
                {
                    y_neg_final[v] += y_neg[v][u];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
            }

            // Handle remaining negative values
            for (; k_neg_loop < end_neg; k_neg_loop++)
            {
                int row_idx = row_index_neg[k_neg_loop];
                for (int v = 0; v < M_UNROLL_FACTOR; v++)
                {
                    y_neg_final[v] += X[(m + v) * K + row_idx];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
            }

            // Final computation and store
            for (int v = 0; v < M_UNROLL_FACTOR; v++)
            {
                Y[(m + v) * N + n] = (y_pos_final[v] - y_neg_final[v]) + b[n];
#ifdef INSTRUMENTATION_RUN
                flops += 2;
#endif
            }
        }
    }

    // Handle remaining M iterations (cleanup loop)
    for (; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            T y_pos[K_UNROLL_FACTOR] = {0};
            T y_neg[K_UNROLL_FACTOR] = {0};

            int k_pos_loop = col_start_pos[n];
            const int end_pos = col_start_pos[n + 1];

            for (; k_pos_loop + K_UNROLL_FACTOR <= end_pos; k_pos_loop += K_UNROLL_FACTOR)
            {
                for (int u = 0; u < K_UNROLL_FACTOR; u++)
                {
                    y_pos[u] += X[m * K + row_index_pos[k_pos_loop + u]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
            }

            T y_pos_final = 0;
            for (int u = 0; u < K_UNROLL_FACTOR; u++)
            {
                y_pos_final += y_pos[u];
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            for (; k_pos_loop < end_pos; k_pos_loop++)
            {
                y_pos_final += X[m * K + row_index_pos[k_pos_loop]];
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            int k_neg_loop = col_start_neg[n];
            const int end_neg = col_start_neg[n + 1];

            for (; k_neg_loop + K_UNROLL_FACTOR <= end_neg; k_neg_loop += K_UNROLL_FACTOR)
            {
                for (int u = 0; u < K_UNROLL_FACTOR; u++)
                {
                    y_neg[u] += X[m * K + row_index_neg[k_neg_loop + u]];
#ifdef INSTRUMENTATION_RUN
                    flops++;
#endif
                }
            }

            T y_neg_final = 0;
            for (int u = 0; u < K_UNROLL_FACTOR; u++)
            {
                y_neg_final += y_neg[u];
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            for (; k_neg_loop < end_neg; k_neg_loop++)
            {
                y_neg_final += X[m * K + row_index_neg[k_neg_loop]];
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            Y[m * N + n] = (y_pos_final - y_neg_final) + b[n];
#ifdef INSTRUMENTATION_RUN
            flops += 2;
#endif
        }
    }
}

#endif