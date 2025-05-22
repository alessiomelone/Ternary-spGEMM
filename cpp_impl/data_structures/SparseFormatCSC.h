#pragma once
#include <vector>
#include "DataStructureInterface.hpp"

/**
 * @brief Compressed Sparse Column (CSC) format for ternary matrices
 *
 * This class implements a CSC format specifically for ternary matrices (-1, 0, 1).
 * It stores positive and negative values separately for better performance.
 */
class SparseFormatCSC : public DataStructureInterface
{
public:
    std::vector<int> col_start_pos;
    std::vector<int> col_start_neg;
    std::vector<int> row_index_pos;
    std::vector<int> row_index_neg;
    int num_rows;
    int num_cols;

    SparseFormatCSC() : num_rows(0), num_cols(0) {}

    SparseFormatCSC(const int *matrix, int K, int N) : num_rows(K), num_cols(N)
    {
        init(matrix, K, N);
    }

    void init(const int *matrix, int rows, int cols) override
    {
        num_rows = rows;
        num_cols = cols;

        // Pre-allocate offset arrays
        col_start_pos.resize(cols + 1, 0);
        col_start_neg.resize(cols + 1, 0);

        // First pass: count non-zeros per column
        for (int n = 0; n < cols; ++n)
        {
            for (int k = 0; k < rows; ++k)
            {
                int val = matrix[k * cols + n];
                if (val == 1)
                {
                    col_start_pos[n + 1]++;
                }
                else if (val == -1)
                {
                    col_start_neg[n + 1]++;
                }
            }
        }

        // Convert counts to start pointers (prefix sum)
        for (int n = 0; n < cols; ++n)
        {
            col_start_pos[n + 1] += col_start_pos[n];
            col_start_neg[n + 1] += col_start_neg[n];
        }

        // Resize index vectors based on total counts
        row_index_pos.resize(col_start_pos[cols]);
        row_index_neg.resize(col_start_neg[cols]);

        // Second pass: fill row indices
        std::vector<int> current_pos_ptr = col_start_pos;
        std::vector<int> current_neg_ptr = col_start_neg;

        for (int n = 0; n < cols; ++n)
        {
            for (int k = 0; k < rows; ++k)
            {
                int val = matrix[k * cols + n];
                if (val == 1)
                {
                    row_index_pos[current_pos_ptr[n]++] = k;
                }
                else if (val == -1)
                {
                    row_index_neg[current_neg_ptr[n]++] = k;
                }
            }
        }
    }

    std::vector<int> getVectorRepresentation(size_t expected_rows, size_t expected_cols) override
    {
        std::vector<int> result(expected_rows * expected_cols, 0);

        for (int n = 0; n < num_cols; ++n)
        {
            // Process positive values
            for (int k = col_start_pos[n]; k < col_start_pos[n + 1]; ++k)
            {
                int row = row_index_pos[k];
                if (row < expected_rows && n < expected_cols)
                {
                    result[row * expected_cols + n] = 1;
                }
            }

            // Process negative values
            for (int k = col_start_neg[n]; k < col_start_neg[n + 1]; ++k)
            {
                int row = row_index_neg[k];
                if (row < expected_rows && n < expected_cols)
                {
                    result[row * expected_cols + n] = -1;
                }
            }
        }

        return result;
    }

    int getNumRows() const override { return num_rows; }
    int getNumCols() const override { return num_cols; }

    void printVars() override
    {
        std::cout << "\nSparseFormatCSC (" << num_rows << "x" << num_cols << "):" << std::endl;
        std::cout << "col_start_pos: ";
        for (int val : col_start_pos)
            std::cout << val << " ";
        std::cout << "\ncol_start_neg: ";
        for (int val : col_start_neg)
            std::cout << val << " ";
        std::cout << "\nrow_index_pos: ";
        for (int val : row_index_pos)
            std::cout << val << " ";
        std::cout << "\nrow_index_neg: ";
        for (int val : row_index_neg)
            std::cout << val << " ";
        std::cout << std::endl;
    }
};