#include "common.h"
#include "SparseGEMM.h" // For SparseFormat
#include "data_structures/CompressedCSC.h"
#include "data_structures/TCSRMatrix.h" // For TCSRMatrix definition

// Rename and modify sparseGEMM_base to be a specific implementation for SparseFormat
template <typename T>
void CSC_base(T *X, const SparseFormat &W_csr, T *b, T *Y, int M, int N, int K)
{
    const int *col_start_pos = W_csr.col_start_pos.data();
    const int *col_start_neg = W_csr.col_start_neg.data();
    const int *row_index_pos = W_csr.row_index_pos.data();
    const int *row_index_neg = W_csr.row_index_neg.data();

    for (int m = 0; m < M; m++)
    {
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

template <typename T>
void CCSC_base(T *X, const CompressedCSC &W, T *b, T *Y, int M, int N, int K)
// X: M rows, K cols
// W: K rows, N cols
// Y: M rows, N cols
{
    const int *col_start = W.col_start.data();
    const int *row_index = W.row_index.data();
    const uint8_t *vals = W.vals.data();

    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            T y_val0 = 0;
            T y_val1 = 0;
            T y_val2 = 0;
            T y_val3 = 0;
            T y_val4 = 0;
            for (int k = col_start[n]; k < col_start[n + 1]; k++)
            {
                const int row = row_index[k];
                const int8_t block = vals[k];
                const int8_t *vals_decoded = decodeCCSC[block];
                y_val0 += vals_decoded[0] * X[m * K + row_index[k] + 0];
                y_val1 += vals_decoded[1] * X[m * K + row_index[k] + 1];
                y_val2 += vals_decoded[2] * X[m * K + row_index[k] + 2];
                y_val3 += vals_decoded[3] * X[m * K + row_index[k] + 3];
                y_val4 += vals_decoded[4] * X[m * K + row_index[k] + 4];
            }
            
            Y[m * N + n + 0] = y_val0 + b[n];
            Y[m * N + n + 1] = y_val1 + b[n];
            Y[m * N + n + 2] = y_val2 + b[n];
            Y[m * N + n + 3] = y_val3 + b[n];
            Y[m * N + n + 4] = y_val4 + b[n];
        }
    }
}

template <typename T>
void TCSR_base(T *X_arg, const TCSRMatrix &W_tcsr, T *B_arg, T *Y_arg,
               int M_dim, int N_dim, int K_dim)
{
    // dimensions check
    // if (W_tcsr.num_matrix_rows != K_dim || W_tcsr.num_matrix_cols != N_dim)
    // {
    //     std::cerr << "TCSR_base Error: Dimension mismatch for W_tcsr." << std::endl;
    //     return;
    // }

    // Initialize Y_arg to 0 (accumulator)
    for (int i = 0; i < M_dim * N_dim; ++i)
    {
        Y_arg[i] = T(0);
    }

    for (int m = 0; m < M_dim; ++m) // Iterate over rows of X and Y
    {
        for (int k = 0; k < K_dim; ++k) // Iterate over columns of X (which are rows of W)
        {
            T x_mk_val = X_arg[m * K_dim + k];
            if (x_mk_val == T(0))
                continue; // Optimization: if X_mk is zero, it contributes nothing

            // Get non-zero elements in row k of W_tcsr
            int row_start_offset_W = W_tcsr.row_offsets[k];
            int row_end_offset_W = W_tcsr.row_offsets[k + 1];

            for (int nz_idx = row_start_offset_W; nz_idx < row_end_offset_W; ++nz_idx)
            {
                std::pair<int, int> decoded_W_element = TCSRMatrix::tcsr_decode_col(W_tcsr.encoded_cols[nz_idx]);
                int n_col_in_W = decoded_W_element.first; // This is the 'n' in W_kn'
                T w_kn_val = static_cast<T>(decoded_W_element.second);

                // Y[m][n_col_in_W] += x_mk_val * w_kn_val;
                Y_arg[m * N_dim + n_col_in_W] += x_mk_val * w_kn_val;
            }
        }
    }

    // Add bias B after Y = XW
    for (int m = 0; m < M_dim; ++m)
    {
        for (int n = 0; n < N_dim; ++n)
        {
            Y_arg[m * N_dim + n] += B_arg[n];
        }
    }
}

// Rename and modify sparseGEMM_unrolled to be a specific implementation for SparseFormat
template <typename T, int UNROLL_FACTOR>
void CSC_unrolled(
    T *X, const SparseFormat &W_csr, T *b, T *Y,
    int M, int N, int K)
{
    const int *col_start_pos = W_csr.col_start_pos.data();
    const int *col_start_neg = W_csr.col_start_neg.data();
    const int *row_index_pos = W_csr.row_index_pos.data();
    const int *row_index_neg = W_csr.row_index_neg.data();

    for (int m = 0; m < M; m++)
    {
        for (int n_idx = 0; n_idx < N; n_idx++)
        { // Renamed n
            T y_pos[UNROLL_FACTOR] = {0};
            T y_neg[UNROLL_FACTOR] = {0};

            int k_pos_loop = col_start_pos[n_idx]; // Renamed k_pos
            const int end_pos = col_start_pos[n_idx + 1];

            for (; k_pos_loop + UNROLL_FACTOR <= end_pos; k_pos_loop += UNROLL_FACTOR)
            {
#pragma unroll
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
                    y_pos[u] += X[m * K + row_index_pos[k_pos_loop + u]];
                }
            }

            T y_pos_final = 0;
            for (int u = 0; u < UNROLL_FACTOR; u++)
            {
                y_pos_final += y_pos[u];
            }

            // remainder loop
            for (; k_pos_loop < end_pos; k_pos_loop++)
            {
                y_pos_final += X[m * K + row_index_pos[k_pos_loop]];
            }

            int k_neg_loop = col_start_neg[n_idx]; // Renamed k_neg
            const int end_neg = col_start_neg[n_idx + 1];

            for (; k_neg_loop + UNROLL_FACTOR <= end_neg; k_neg_loop += UNROLL_FACTOR)
            {
#pragma unroll
                for (int u = 0; u < UNROLL_FACTOR; u++)
                {
                    y_neg[u] += X[m * K + row_index_neg[k_neg_loop + u]];
                }
            }

            T y_neg_final = 0;
            for (int u = 0; u < UNROLL_FACTOR; u++)
            {
                y_neg_final += y_neg[u];
            }

            // remainder loop
            for (; k_neg_loop < end_neg; k_neg_loop++)
            {
                y_neg_final += X[m * K + row_index_neg[k_neg_loop]];
            }

            Y[m * N + n_idx] = (y_pos_final - y_neg_final) + b[n_idx];
        }
    }
}

// --- Explicit Instantiations ---
// This tells the compiler to generate code for these specific versions in comp.o
template void CSC_base<float>(float *, const SparseFormat &, float *, float *, int, int, int);
template void CCSC_base<float>(float *, const CompressedCSC &, float *, float *, int, int, int);
template void TCSR_base<float>(float *, const TCSRMatrix &, float *, float *, int, int, int);
template void CSC_unrolled<float, 2>(float *, const SparseFormat &, float *, float *, int, int, int);
// If you use other unroll factors or other types for T, you'd add them here.
template void CSC_unrolled<float, 12>(float *, const SparseFormat &, float *, float *, int, int, int);
