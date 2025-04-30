// Compile and run: $ gcc -O3 ref.c && ./a.out
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct ternarySparseFormat_t {
	int *col_start_pos;
	int *col_start_neg;
	int *row_index_pos;
	int *row_index_neg;
	int size; // TODO: Check if we the size of these arrays can be reduced from N
} ternarySparseFormat;

// TODO: Run valgrind on it to verify no unsafe mem accesses -- should be ok since tests pass
bool compare_results(int* result, int* groundTruth, int H, int W) {
	for (int h = 0; h < H; h++) {
		for (int w = 0; w < W; w++) {
			int i = h * W + w;
			if (abs(result[i] - groundTruth[i]) > 10e-6) {
				printf("Error at: H=%d, W=%d, result=%d, groundTruth=%d\n", h, w, result[i], groundTruth[i]);
				return false;
			}
		}
	}

	return true;
}

// Do basic GEMM (3-fold loop)
void GEMM(int* X, int* W, int* b, int* Y, int M, int N, int K) {
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			int y = 0;
			for (int k = 0; k < K; k++) {
				y += X[m * K + k] * W[k * N + n];
			}
			Y[m * N + n] = y + b[n];
		}
	}
}

// Do Sparse GEMM, store results in parameter Y
void sparseGEMM(int* X, ternarySparseFormat *W, int* b, int* Y, int M, int N, int K) {
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			int y = 0;
			for (int k = W->col_start_pos[n]; k < W->col_start_pos[n + 1]; k++) {
				y += X[m * K + W->row_index_pos[k]];
			}
			for (int k = W->col_start_neg[n]; k < W->col_start_neg[n + 1]; k++) {
				y -= X[m * K + W->row_index_neg[k]];
			}
			Y[m * N + n] = y + b[n];
		}
	}
}

// Convert ternary matrix to Ternary Sparce Format
ternarySparseFormat *convertTernaryToSparseFormat(int* matrix, int K, int N, int nonZeroPercentage) {
    // TODO: Verify sizes of each sub array
    int nonZeroVals = (K * N) / (double) (nonZeroPercentage) + 1;
	ternarySparseFormat *tsf = malloc(sizeof(ternarySparseFormat));
	tsf->size = N+1;
	tsf->col_start_pos = malloc(tsf->size * sizeof(int));
	tsf->col_start_neg = malloc(tsf->size * sizeof(int));
	tsf->row_index_pos = malloc(nonZeroVals * sizeof(int));
	tsf->row_index_neg = malloc(nonZeroVals * sizeof(int));

	int column_start_pos = 0;
	int column_start_neg = 0;
	int row_index_pos_ind = 0;
	int row_index_neg_ind = 0;
	int n;
	for (n = 0; n < N; n++) {
		tsf->col_start_pos[n] = column_start_pos;
		tsf->col_start_neg[n] = column_start_neg;
		for (int k = 0; k < K; k++) {
			if (matrix[k * N + n] >= 1) {
				column_start_pos++;
				tsf->row_index_pos[row_index_pos_ind++] = k;
			}
			else if (matrix[k * N + n] <= -1) {
				column_start_neg++;
				tsf->row_index_neg[row_index_neg_ind++] = k;
			}
		}
	}
	tsf->col_start_pos[n] = column_start_pos;
	tsf->col_start_neg[n] = column_start_neg;
    return tsf;
}

// Free memory from malloc()
void destroyTernarySparceFormat(ternarySparseFormat *tsf) {
	free(tsf->col_start_pos);
	free(tsf->col_start_neg);
	free(tsf->row_index_pos);
	free(tsf->row_index_neg);
	free(tsf);
}


//// TEST INFRA


int *generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution) {
    // TODO : Free y
    int *y = malloc(sizeof(int) * H * W);
    if (uniformDistribution) {
        for (int h = 0; h < H; h++) {
            for (int w = 0; w < W; w += nonZero * 2) {
                // Assign +1, -1 to each 2 x nonZero slots
                int randomA = rand() % nonZero * 2;
                int randomB = rand() % nonZero * 2;
                y[w + randomA] = 1;
                while (randomA==randomB) {
                    randomB = rand() % nonZero * 2;
                }
                y[w + randomB] = 1;
            }
        }
    }
    else {
        for (int h = 0; h < H; h++) {
            // Assign +1 to W / nonZero / 2 places
            int count = 0;
            int limit = (W / nonZero) / 2 - 1;
            while (count < limit) {
                int randomA = rand() % W;
                if (y[randomA] == 0) {
                    y[randomA] = 1;
                    count++;
                }
            }

            // Assign -1 to W / nonZero / 2 places
            count = 0;
            while (count < (W / nonZero / 2 - 1)) {
                int randomA = rand() % W;
                if (y[randomA] == 0) {
                    y[randomA] = -1;
                    count++;
                }
            }
        }
    }

    return y;
}

int main() {
    int TEST_CASES = 8;
    int M[] = { 1, 16, 64, 256, 1000, 4000, 16000, 64000 };
    int K[] = { 512, 1024, 2048,  4096, 2048, 4096, 8192, 16384 };
    int N[] = { 2048, 4096, 8192, 16384,  512, 1024, 2048,  4096 };
    int nonZero[] = { 2, 4, 8, 16 }; // 1/2, 1/4, 1/8, 1/16 non-zero values

    int * X = generateSparseMatrix(M[7], K[7], 4, true);
    // cout << "X initialized."<< endl;
    int * W = generateSparseMatrix(K[7], N[7], nonZero[3], true);
    // cout << "W initialized." << endl;

    for (int m = 0; m < TEST_CASES; m++) {
        for (int n = 0; n < TEST_CASES; n++) {
            printf("M=%d, N=%d, K=%d\n", M[m], N[n], K[n]);
            int *Y = calloc(M[m] * N[n], sizeof(int));
            int *B = malloc(N[n] * sizeof(int));
            for (int i = 0; i < N[n]; i++) {
                B[i] = 2;
            }
            int *refY = calloc(M[m] * N[n], sizeof(int));
            ternarySparseFormat *tsf = convertTernaryToSparseFormat(W, K[n], N[n], nonZero[3]);
            sparseGEMM(X, tsf, B, Y, M[m], N[n], K[n]);
            GEMM(X, W, B, refY, M[m], N[n], K[n]);
            if (compare_results(Y, refY, M[m], N[n])) {
                printf("Test case %d out of %d passed!\n", m * TEST_CASES + n + 1, TEST_CASES * TEST_CASES);
            }
            destroyTernarySparceFormat(tsf);
            free(refY);
            free(B);
            free(Y);
        }
    }

    free(X);
    free(W);
}