#include "common.h"
#include "SparseGEMM.h" // For SparseFormat
#include "data_structures/CompressedCSC.h"
#include "data_structures/TCSRMatrix.h" // For TCSRMatrix definition
#include "data_structures/TCSCMatrix.h"

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
    const int *row_offsets_data = W_tcsr.row_offsets.data();
    const int *encoded_cols_data = W_tcsr.encoded_cols.data();

    for (int m = 0; m < M_dim; ++m)
    {
        const T *X_row_m = X_arg + m * K_dim;
        T *Y_row_m = Y_arg + m * N_dim;

        // Initialize current row of Y with bias
        for (int n = 0; n < N_dim; ++n)
        {
            Y_row_m[n] = B_arg[n];
        }

        for (int k = 0; k < K_dim; ++k)
        {
            T x_mk_val = X_row_m[k];

            int row_start_offset_W = row_offsets_data[k];
            int row_end_offset_W = row_offsets_data[k + 1];

            for (int nz_idx = row_start_offset_W; nz_idx < row_end_offset_W; ++nz_idx)
            {
                std::pair<int, int> decoded_W_element = TCSRMatrix::tcsr_decode_col(encoded_cols_data[nz_idx]);
                int n_col_in_W = decoded_W_element.first;
                T w_kn_val = static_cast<T>(decoded_W_element.second);

                Y_row_m[n_col_in_W] += x_mk_val * w_kn_val;
            }
        }
    }
}

template <typename T, int UNROLL_FACTOR>
void TCSR_unrolled(T *X_arg, const TCSRMatrix &W_tcsr, T *B_arg, T *Y_arg,
                   int M_dim, int N_dim, int K_dim)
{
    const int *row_offsets_data = W_tcsr.row_offsets.data();
    const int *encoded_cols_data = W_tcsr.encoded_cols.data();

    for (int m = 0; m < M_dim; ++m)
    {
        const T *X_row_m = X_arg + m * K_dim;
        T *Y_row_m = Y_arg + m * N_dim;

        // Initialize current row of Y with bias
        for (int n_init = 0; n_init < N_dim; ++n_init)
        {
            Y_row_m[n_init] = B_arg[n_init];
        }

        for (int k = 0; k < K_dim; ++k)
        {
            T x_mk_val = X_row_m[k];

            const int row_start_offset_W = row_offsets_data[k];
            const int row_end_offset_W = row_offsets_data[k + 1];
            int nz_idx = row_start_offset_W;

            // Unrolled part
            for (; nz_idx + UNROLL_FACTOR <= row_end_offset_W; nz_idx += UNROLL_FACTOR)
            {
#pragma unroll
                for (int u = 0; u < UNROLL_FACTOR; ++u)
                {
                    std::pair<int, int> decoded_W_element = TCSRMatrix::tcsr_decode_col(encoded_cols_data[nz_idx + u]);
                    int n_col_in_W = decoded_W_element.first;
                    T w_kn_val = static_cast<T>(decoded_W_element.second);
                    Y_row_m[n_col_in_W] += x_mk_val * w_kn_val;
                }
            }

            // Remainder loop
            for (; nz_idx < row_end_offset_W; ++nz_idx)
            {
                std::pair<int, int> decoded_W_element = TCSRMatrix::tcsr_decode_col(encoded_cols_data[nz_idx]);
                int n_col_in_W = decoded_W_element.first;
                T w_kn_val = static_cast<T>(decoded_W_element.second);
                Y_row_m[n_col_in_W] += x_mk_val * w_kn_val;
            }
        }
    }
}

template <typename T>
void TCSC_base(T *X_arg, const TCSCMatrix &W_tcsc, T *B_arg, T *Y_arg,
               int M_dim, int N_dim, int K_dim)
{
    const int *col_offsets_data = W_tcsc.col_offsets.data();
    const int *encoded_rows_data = W_tcsc.encoded_rows.data();

    for (int m = 0; m < M_dim; ++m)
    {
        const T *X_row_m = X_arg + m * K_dim;
        T *Y_row_m = Y_arg + m * N_dim;

        for (int n = 0; n < N_dim; ++n)
        {
            T y_mn_acc = 0;
            int col_start_offset_W = col_offsets_data[n];
            int col_end_offset_W = col_offsets_data[n + 1];

            for (int nz_idx = col_start_offset_W; nz_idx < col_end_offset_W; ++nz_idx)
            {
                std::pair<int, int> decoded_W = TCSCMatrix::tcsc_decode_row(encoded_rows_data[nz_idx]);
                int k_row_in_W = decoded_W.first;
                T w_kn_val = static_cast<T>(decoded_W.second);
                y_mn_acc += X_row_m[k_row_in_W] * w_kn_val;
            }
            Y_row_m[n] = y_mn_acc + B_arg[n];
        }
    }
}

template <typename T, int UNROLL_FACTOR>
void TCSC_unrolled(T *X_arg, const TCSCMatrix &W_tcsc, T *B_arg, T *Y_arg,
                   int M_dim, int N_dim, int K_dim)
{
    const int *encoded_rows_data = W_tcsc.encoded_rows.data();
    const int *col_offsets_data = W_tcsc.col_offsets.data();

    // For each row m of X and Y
    for (int m = 0; m < M_dim; ++m)
    {
        const T *X_row_m = X_arg + m * K_dim; // Pointer to current row of X
        T *Y_row_m = Y_arg + m * N_dim;       // Pointer to current row of Y

        // For each column n of W and Y
        for (int n = 0; n < N_dim; ++n)
        {
            T y_mn_acc_partials[UNROLL_FACTOR];
            // Initialize partial accumulators to zero
            for (int u = 0; u < UNROLL_FACTOR; ++u)
            {
                y_mn_acc_partials[u] = T(0);
            }

            const int col_start_offset_W = col_offsets_data[n];
            const int col_end_offset_W = col_offsets_data[n + 1];
            int nz_idx = col_start_offset_W;

            for (; nz_idx + UNROLL_FACTOR <= col_end_offset_W; nz_idx += UNROLL_FACTOR)
            {

#pragma unroll
                for (int u = 0; u < UNROLL_FACTOR; ++u)
                {
                    // Decode W element (row index in W and value)
                    std::pair<int, int> decoded_W = TCSCMatrix::tcsc_decode_row(encoded_rows_data[nz_idx + u]);
                    int k_row_in_W = decoded_W.first;              // This is the 'k' index for X
                    T w_kn_val = static_cast<T>(decoded_W.second); // This is the W[k,n] value

                    // Accumulate product into the corresponding partial accumulator
                    y_mn_acc_partials[u] += X_row_m[k_row_in_W] * w_kn_val;
                }
            }

            // Sum up partial accumulators
            T y_mn_acc_final = T(0);
            for (int u = 0; u < UNROLL_FACTOR; ++u)
            {
                y_mn_acc_final += y_mn_acc_partials[u];
            }

            // Remainder loop
            for (; nz_idx < col_end_offset_W; ++nz_idx)
            {
                std::pair<int, int> decoded_W = TCSCMatrix::tcsc_decode_row(encoded_rows_data[nz_idx]);
                int k_row_in_W = decoded_W.first;
                T w_kn_val = static_cast<T>(decoded_W.second);
                y_mn_acc_final += X_row_m[k_row_in_W] * w_kn_val;
            }

            Y_row_m[n] = y_mn_acc_final + B_arg[n];
        }
    }
}

template <typename T, int UNROLL_FACTOR, int TILE_M, int TILE_N>
void TCSC_unrolled_tiled(T *X_arg, const TCSCMatrix &W_tcsc, T *B_arg, T *Y_arg,
                         int M_dim, int N_dim, int K_dim)
{
    const int *encoded_rows_data = W_tcsc.encoded_rows.data();
    const int *col_offsets_data = W_tcsc.col_offsets.data();

    for (int m_tile_start = 0; m_tile_start < M_dim; m_tile_start += TILE_M)
    {
        int m_tile_end = std::min(m_tile_start + TILE_M, M_dim);

        for (int n_tile_start = 0; n_tile_start < N_dim; n_tile_start += TILE_N)
        {
            int n_tile_end = std::min(n_tile_start + TILE_N, N_dim);

            // Process tile: Y[m_tile_start:m_tile_end, n_tile_start:n_tile_end]
            for (int m = m_tile_start; m < m_tile_end; ++m)
            {
                const T *X_row_m = X_arg + m * K_dim;
                T *Y_row_m = Y_arg + m * N_dim;

                for (int n = n_tile_start; n < n_tile_end; ++n)
                {
                    T y_mn_acc_partials[UNROLL_FACTOR];
                    for (int u = 0; u < UNROLL_FACTOR; ++u)
                    {
                        y_mn_acc_partials[u] = T(0);
                    }

                    const int col_start_offset_W = col_offsets_data[n];
                    const int col_end_offset_W = col_offsets_data[n + 1];
                    int nz_idx = col_start_offset_W;

                    // Unrolled part
                    for (; nz_idx + UNROLL_FACTOR <= col_end_offset_W; nz_idx += UNROLL_FACTOR)
                    {
#pragma unroll
                        for (int u = 0; u < UNROLL_FACTOR; ++u)
                        {
                            std::pair<int, int> decoded_W = TCSCMatrix::tcsc_decode_row(encoded_rows_data[nz_idx + u]);
                            int k_row_in_W = decoded_W.first;
                            T w_kn_val = static_cast<T>(decoded_W.second);
                            y_mn_acc_partials[u] += X_row_m[k_row_in_W] * w_kn_val;
                        }
                    }

                    T y_mn_acc_final = T(0);
                    for (int u = 0; u < UNROLL_FACTOR; ++u)
                    {
                        y_mn_acc_final += y_mn_acc_partials[u];
                    }

                    // Remainder loop
                    for (; nz_idx < col_end_offset_W; ++nz_idx)
                    {
                        std::pair<int, int> decoded_W = TCSCMatrix::tcsc_decode_row(encoded_rows_data[nz_idx]);
                        int k_row_in_W = decoded_W.first;
                        T w_kn_val = static_cast<T>(decoded_W.second);
                        y_mn_acc_final += X_row_m[k_row_in_W] * w_kn_val;
                    }
                    Y_row_m[n] = y_mn_acc_final + B_arg[n];
                }
            }
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
template void TCSC_base<float>(float *, const TCSCMatrix &, float *, float *, int, int, int);
template void CSC_unrolled<float, 2>(float *, const SparseFormat &, float *, float *, int, int, int);
// If you use other unroll factors or other types for T, you'd add them here.
template void CSC_unrolled<float, 12>(float *, const SparseFormat &, float *, float *, int, int, int);

template void TCSR_unrolled<float, 12>(float *, const TCSRMatrix &, float *, float *, int, int, int);
template void TCSC_unrolled<float, 12>(float *, const TCSCMatrix &, float *, float *, int, int, int);

template void TCSC_unrolled_tiled<float, 12, 32, 32>(float *, const TCSCMatrix &, float *, float *, int, int, int);