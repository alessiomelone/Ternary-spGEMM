#ifndef COMP_H_PRELU
#define COMP_H_PRELU

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

template <typename T>
void BaseTCSC_PreLU(T *X, const TCSC &W_csc, T *b, T *alpha, T *Y, int M, int N, int K)
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

            // Apply bias
            y_val = y_val + b[n];
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif

            // Apply PReLU activation: f(x) = x if x > 0, f(x) = α * x if x ≤ 0
            if (y_val > 0)
            {
                Y[m * N + n] = y_val;
            }
            else
            {
                Y[m * N + n] = alpha[n] * y_val;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }
        }
    }
}

#endif // COMP_H_PRELU
