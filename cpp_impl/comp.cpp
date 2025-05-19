#include "common.h"
#include "SparseGEMM.h" // For SparseFormat

// Rename and modify sparseGEMM_base to be a specific implementation for SparseFormat
template <typename T>
void sparseGEMM_csr_base_impl(T *X, const SparseFormat& W_csr, T *b, T *Y, int M, int N, int K)
{
    const int* col_start_pos = W_csr.col_start_pos.data();
    const int* col_start_neg = W_csr.col_start_neg.data();
    const int* row_index_pos = W_csr.row_index_pos.data();
    const int* row_index_neg = W_csr.row_index_neg.data();

    for (int m = 0; m < M; m++)
    {
        for (int n_idx = 0; n_idx < N; n_idx++) // Renamed n to n_idx
        {
            T y_val = 0; // Renamed y to y_val
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

// Rename and modify sparseGEMM_unrolled to be a specific implementation for SparseFormat
template <typename T, int UNROLL_FACTOR>
void sparseGEMM_csr_unrolled_impl(
    T *X, const SparseFormat& W_csr, T *b, T *Y,
    int M, int N, int K)
{
    const int* col_start_pos = W_csr.col_start_pos.data();
    const int* col_start_neg = W_csr.col_start_neg.data();
    const int* row_index_pos = W_csr.row_index_pos.data();
    const int* row_index_neg = W_csr.row_index_neg.data();

    for (int m = 0; m < M; m++) {
        for (int n_idx = 0; n_idx < N; n_idx++) { // Renamed n
            T y_pos[UNROLL_FACTOR] = {0};
            T y_neg[UNROLL_FACTOR] = {0};

            int k_pos_loop = col_start_pos[n_idx]; // Renamed k_pos
            const int end_pos = col_start_pos[n_idx + 1];
            
            for (; k_pos_loop + UNROLL_FACTOR <= end_pos; k_pos_loop += UNROLL_FACTOR) {
                #pragma unroll
                for (int u = 0; u < UNROLL_FACTOR; u++) {
                    y_pos[u] += X[m * K + row_index_pos[k_pos_loop + u]];
                }
            }

            T y_pos_final = 0;
            for (int u = 0; u < UNROLL_FACTOR; u++) {
                y_pos_final += y_pos[u];
            }

            // remainder loop
            for (; k_pos_loop < end_pos; k_pos_loop++) {
                y_pos_final += X[m * K + row_index_pos[k_pos_loop]];
            }

            int k_neg_loop = col_start_neg[n_idx]; // Renamed k_neg
            const int end_neg = col_start_neg[n_idx + 1];
            
            for (; k_neg_loop + UNROLL_FACTOR <= end_neg; k_neg_loop += UNROLL_FACTOR) {
                #pragma unroll
                for (int u = 0; u < UNROLL_FACTOR; u++) {
                    y_neg[u] += X[m * K + row_index_neg[k_neg_loop + u]];
                }
            }

            T y_neg_final = 0;
            for (int u = 0; u < UNROLL_FACTOR; u++) {
                y_neg_final += y_neg[u];
            }

            // remainder loop
            for (; k_neg_loop < end_neg; k_neg_loop++) {
                y_neg_final += X[m * K + row_index_neg[k_neg_loop]];
            }

            Y[m * N + n_idx] = (y_pos_final - y_neg_final) + b[n_idx];
        }
    }
}


// --- Explicit Instantiations ---
// This tells the compiler to generate code for these specific versions in comp.o
template void sparseGEMM_csr_base_impl<float>(float*, const SparseFormat&, float*, float*, int, int, int);
template void sparseGEMM_csr_unrolled_impl<float, 2>(float*, const SparseFormat&, float*, float*, int, int, int);
// If you use other unroll factors or other types for T, you'd add them here.
template void sparseGEMM_csr_unrolled_impl<float, 12>(float*, const SparseFormat&, float*, float*, int, int, int);