#include "sparse_format.h"

// Do Sparse GEMM, store results in parameter Y
void sparseGEMM(int *X, ternarySparseFormat *W, int *b, int *Y, int M, int N, int K)
{
    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            int y = 0;
            for (int k = W->col_start_pos[n]; k < W->col_start_pos[n + 1]; k++)
            {
                y += X[m * K + W->row_index_pos[k]];
            }
            for (int k = W->col_start_neg[n]; k < W->col_start_neg[n + 1]; k++)
            {
                y -= X[m * K + W->row_index_neg[k]];
            }
            Y[m * N + n] = y + b[n];
        }
    }
}

void GEMM(int *X, int *W, int *b, int *Y, int M, int N, int K)
{
    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            int y = 0;
            for (int k = 0; k < K; k++)
            {
                y += X[m * K + k] * W[k * N + n];
            }
            Y[m * N + n] = y + b[n];
        }
    }
}

// Convert ternary matrix to Ternary Sparse Format
ternarySparseFormat *convertTernaryToSparseFormat(int *matrix, int K, int N, int nonZeroPercentage)
{
    // TODO: Verify sizes of each sub array
    int nonZeroVals = (K * N) / (double)(nonZeroPercentage) + 1;
    ternarySparseFormat *tsf = malloc(sizeof(ternarySparseFormat));
    tsf->size = N + 1;
    tsf->col_start_pos = malloc(tsf->size * sizeof(int));
    tsf->col_start_neg = malloc(tsf->size * sizeof(int));
    tsf->row_index_pos = malloc(nonZeroVals * sizeof(int));
    tsf->row_index_neg = malloc(nonZeroVals * sizeof(int));

    int column_start_pos = 0;
    int column_start_neg = 0;
    int row_index_pos_ind = 0;
    int row_index_neg_ind = 0;
    int n;
    for (n = 0; n < N; n++)
    {
        tsf->col_start_pos[n] = column_start_pos;
        tsf->col_start_neg[n] = column_start_neg;
        for (int k = 0; k < K; k++)
        {
            if (matrix[k * N + n] >= 1)
            {
                column_start_pos++;
                tsf->row_index_pos[row_index_pos_ind++] = k;
            }
            else if (matrix[k * N + n] <= -1)
            {
                column_start_neg++;
                tsf->row_index_neg[row_index_neg_ind++] = k;
            }
        }
    }
    tsf->col_start_pos[n] = column_start_pos;
    tsf->col_start_neg[n] = column_start_neg;
    return tsf;
}

// TODO: Run valgrind on it to verify no unsafe mem accesses -- should be ok since tests pass
bool compare_results(int *result, int *groundTruth, int H, int W)
{
    for (int h = 0; h < H; h++)
    {
        for (int w = 0; w < W; w++)
        {
            int i = h * W + w;
            if (abs(result[i] - groundTruth[i]) > 10e-6)
            {
                printf("Error at: H=%d, W=%d, result=%d, groundTruth=%d\n", h, w, result[i], groundTruth[i]);
                return false;
            }
        }
    }

    return true;
}


float *generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution)
{
    long long totalElements = (long long)H * W;
    float *y = (float *)malloc(sizeof(float) * totalElements);
    memset(y, 0, sizeof(float) * totalElements);

    // Assuming nonZero > 0
    long long numNonZeroTarget = totalElements / nonZero;

    if (uniformDistribution)
    {
        long long numPos = numNonZeroTarget / 2;
        long long numNeg = numNonZeroTarget - numPos;
        long long count = 0;

        while (count < numPos) {
            long long index = rand() % totalElements;
            if (y[index] == 0.0f) {
                y[index] = 1.0f;
                count++;
            }
        }
        count = 0;
        while (count < numNeg) {
            long long index = rand() % totalElements;
            if (y[index] == 0.0f) {
                y[index] = -1.0f;
                count++;
            }
        }
    }
    else
    {
        for (int h = 0; h < H; ++h) {
            int numNonZeroForRow = W / nonZero;
             if (numNonZeroForRow > W) numNonZeroForRow = W;

            int numPosForRow = numNonZeroForRow / 2;
            int numNegForRow = numNonZeroForRow - numPosForRow;
            int count = 0;

            while (count < numPosForRow) {
                int randomCol = rand() % W;
                long long index = (long long)h * W + randomCol;
                if (y[index] == 0.0f) {
                    y[index] = 1.0f;
                    count++;
                }
            }
            count = 0;
             while (count < numNegForRow) {
                int randomCol = rand() % W;
                 long long index = (long long)h * W + randomCol;
                if (y[index] == 0.0f) {
                    y[index] = -1.0f;
                    count++;
                }
            }
        }
    }
    return y; // Remember caller must free this memory.
}

// Free memory for TSF
void destroyTernarySparcseFormat(ternarySparseFormat *tsf)
{
    free(tsf->col_start_pos);
    free(tsf->col_start_neg);
    free(tsf->row_index_pos);
    free(tsf->row_index_neg);
    free(tsf);
}