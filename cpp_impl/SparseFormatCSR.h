#ifndef SPARSEFORMATCSR_H
#define SPARSEFORMATCSR_H

#include <vector>
#include <numeric> // For std::iota if needed, or explicit loops

// Assumes W_raw is a KxN matrix stored in row-major order.
struct SparseFormatCSR {
    std::vector<int> row_start_pos;
    std::vector<int> row_start_neg;
    std::vector<int> col_index_pos;
    std::vector<int> col_index_neg;
    int num_rows; // K
    int num_cols; // N

    SparseFormatCSR(const int *W_raw, int K_rows, int N_cols) : num_rows(K_rows), num_cols(N_cols) {
        row_start_pos.resize(K_rows + 1, 0);
        row_start_neg.resize(K_rows + 1, 0);

        // First pass: count non-zeros per row to determine sizes for col_index vectors
        // and to populate row_start arrays (as counts initially)
        for (int k = 0; k < K_rows; ++k) {
            for (int n = 0; n < N_cols; ++n) {
                int val = W_raw[k * N_cols + n];
                if (val == 1) {
                    row_start_pos[k + 1]++;
                } else if (val == -1) {
                    row_start_neg[k + 1]++;
                }
            }
        }

        // Convert counts in row_start to actual start pointers (prefix sum)
        for (int k = 0; k < K_rows; ++k) {
            row_start_pos[k + 1] += row_start_pos[k];
            row_start_neg[k + 1] += row_start_neg[k];
        }

        // Resize col_index vectors based on total counts
        col_index_pos.resize(row_start_pos[K_rows]);
        col_index_neg.resize(row_start_neg[K_rows]);

        // Second pass: fill col_index vectors
        // Need temporary current_pos arrays, similar to how CSC does it, or use the row_start arrays carefully
        std::vector<int> current_pos_ptr = row_start_pos; // copy for current insertion points
        std::vector<int> current_neg_ptr = row_start_neg; // copy for current insertion points

        for (int k = 0; k < K_rows; ++k) {
            for (int n = 0; n < N_cols; ++n) {
                int val = W_raw[k * N_cols + n];
                if (val == 1) {
                    col_index_pos[current_pos_ptr[k]++] = n;
                } else if (val == -1) {
                    col_index_neg[current_neg_ptr[k]++] = n;
                }
            }
        }
    }
};

#endif // SPARSEFORMATCSR_H 