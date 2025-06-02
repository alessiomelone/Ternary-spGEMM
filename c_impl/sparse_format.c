#include "sparse_format.h"

// Do Sparse GEMM, store results in parameter Y
void sparseGEMM(float *X, ternarySparseFormat *W, float *b, float *Y, int M, int N, int K)
{
    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            float y = 0;
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

void GEMM(float *X, int *W, float *b, float *Y, int M, int N, int K)
{
    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            float y = 0;
            for (int k = 0; k < K; k++)
            {
                y += X[m * K + k] * W[k * N + n];
            }
            Y[m * N + n] = y + b[n];
        }
    }
}


float *initX(int LEN, int Range) {
	float *X = malloc(LEN * sizeof(float));
	for (int i = 0; i < LEN; i++) {
        float random_value = ((float)rand() / (float)RAND_MAX) * (2.0f * Range) - Range;
		X[i] = random_value;
	}
	return X;
};

// Convert ternary matrix to Ternary Sparse Format
ternarySparseFormat *convertTernaryToSparseFormat(int *matrix, int K, int N, int nonZero)
{
    // TODO: Verify sizes of each sub array
    int nonZeroVals = (K * N) / (float)(nonZero) + 1;
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
            if (matrix[k * N + n] == 1)
            {
                column_start_pos++;
                tsf->row_index_pos[row_index_pos_ind++] = k;
            }
            else if (matrix[k * N + n] == -1)
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
bool compare_results(float *result, float *groundTruth, int H, int W)
{
    for (int h = 0; h < H; h++)
    {
        for (int w = 0; w < W; w++)
        {
            int i = h * W + w;
            if (fabsf(result[i] - groundTruth[i]) > (float) ERROR_TOLERANCE)
            {
                printf("Error at: H=%d, W=%d, result=%f, groundTruth=%f\n", h, w, result[i], groundTruth[i]);
                return false;
            }
        }
    }

    return true;
}


int *generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution) {
    int *y = calloc(H * W, sizeof(int));
    if (!y) return NULL;

    /* Uniform distribution: process each row in blocks of size 2*nonZero */
    if (uniformDistribution) {
        for (int h = 0; h < H; ++h) {
            for (int w = 0; w < W; w += nonZero * 2) {
                int segmentSize = nonZero * 2;
                if (w + segmentSize > W)  /* last partial segment if any */
                    segmentSize = W - w;

                int a = rand() % segmentSize;
                int b = rand() % segmentSize;
                while (b == a)
                    b = rand() % segmentSize;

                y[h * W + w + a] = 1;
                y[h * W + w + b] = -1;
            }
        }
    }
    /* Non-uniform: random +/- per row with a small variation */
    else {
        int basePerRow = W / nonZero;
        int variMax    = basePerRow / 20 + 1;

        for (int h = 0; h < H; ++h) {
            int posVari   = rand() % (variMax + 1);
            int limitPos  = basePerRow / 2 + posVari;
            int limitNeg  = basePerRow / 2 - posVari;
            int count;

            /* place +1 */
            count = 0;
            while (count < limitPos) {
                int idx = rand() % W;
                if (y[h * W + idx] == 0) {
                    y[h * W + idx] = 1;
                    ++count;
                }
            }
            /* place -1 */
            count = 0;
            while (count < limitNeg) {
                int idx = rand() % W;
                if (y[h * W + idx] == 0) {
                    y[h * W + idx] = -1;
                    ++count;
                }
            }
        }
    }

    return y;
}

// Free memory for TSF
void destroyTernarySparseFormat(ternarySparseFormat *tsf)
{
    free(tsf->col_start_pos);
    free(tsf->col_start_neg);
    free(tsf->row_index_pos);
    free(tsf->row_index_neg);
    free(tsf);
}