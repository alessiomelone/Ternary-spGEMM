#pragma once
#include <vector>
#include <iostream>
#include "DataStructureInterface.hpp"

/**
 * @brief Compressed Sparse Row (CSR) format for ternary matrices
 *
 * This class implements a CSR format specifically for ternary matrices (-1, 0, 1).
 * It stores positive and negative values separately for better performance.
 */
class BaseTCSR : public DataStructureInterface
{
public:
    std::vector<int> row_start_pos;
    std::vector<int> row_start_neg;
    std::vector<int> col_index_pos;
    std::vector<int> col_index_neg;
    int num_rows;
    int num_cols;

    BaseTCSR() : num_rows(0), num_cols(0) {}

    BaseTCSR(const int *matrix, int K, int N) : num_rows(K), num_cols(N)
    {
        init(matrix, K, N);
    }

    void init(const int *matrix, int rows, int cols) override
    {
        num_rows = rows;
        num_cols = cols;

        // Pre-allocate offset arrays
        row_start_pos.resize(rows + 1, 0);
        row_start_neg.resize(rows + 1, 0);

        // First pass: count non-zeros per row
        for (int k = 0; k < rows; ++k)
        {
            for (int n = 0; n < cols; ++n)
            {
                int val = matrix[k * cols + n];
                if (val == 1)
                {
                    row_start_pos[k + 1]++;
                }
                else if (val == -1)
                {
                    row_start_neg[k + 1]++;
                }
            }
        }

        // Convert counts to start pointers (prefix sum)
        for (int k = 0; k < rows; ++k)
        {
            row_start_pos[k + 1] += row_start_pos[k];
            row_start_neg[k + 1] += row_start_neg[k];
        }

        // Resize index vectors based on total counts
        col_index_pos.resize(row_start_pos[rows]);
        col_index_neg.resize(row_start_neg[rows]);

        // Second pass: fill column indices
        std::vector<int> current_pos_ptr = row_start_pos;
        std::vector<int> current_neg_ptr = row_start_neg;

        for (int k = 0; k < rows; ++k)
        {
            for (int n = 0; n < cols; ++n)
            {
                int val = matrix[k * cols + n];
                if (val == 1)
                {
                    col_index_pos[current_pos_ptr[k]++] = n;
                }
                else if (val == -1)
                {
                    col_index_neg[current_neg_ptr[k]++] = n;
                }
            }
        }
    }

    std::vector<int> getVectorRepresentation(size_t expected_rows, size_t expected_cols) override
    {
        std::vector<int> result(expected_rows * expected_cols, 0);

        for (int r = 0; r < num_rows; ++r)
        {
            // Process positive values
            for (int k = row_start_pos[r]; k < row_start_pos[r + 1]; ++k)
            {
                int col = col_index_pos[k];
                if (r < expected_rows && col < expected_cols)
                {
                    result[r * expected_cols + col] = 1;
                }
            }

            // Process negative values
            for (int k = row_start_neg[r]; k < row_start_neg[r + 1]; ++k)
            {
                int col = col_index_neg[k];
                if (r < expected_rows && col < expected_cols)
                {
                    result[r * expected_cols + col] = -1;
                }
            }
        }

        return result;
    }

    int getNumRows() const { return num_rows; }
    int getNumCols() const { return num_cols; }

    void printVars() override
    {
        std::cout << "\nSparseFormatCSR (" << num_rows << "x" << num_cols << "):" << std::endl;
        std::cout << "row_start_pos: ";
        for (int val : row_start_pos)
            std::cout << val << " ";
        std::cout << "\nrow_start_neg: ";
        for (int val : row_start_neg)
            std::cout << val << " ";
        std::cout << "\ncol_index_pos: ";
        for (int val : col_index_pos)
            std::cout << val << " ";
        std::cout << "\ncol_index_neg: ";
        for (int val : col_index_neg)
            std::cout << val << " ";
        std::cout << std::endl;
    }

    int getDataStructureSize() const
    {
        return sizeof(int) * (row_start_pos.size() +
                              row_start_neg.size() +
                              col_index_pos.size() +
                              col_index_neg.size() +
                              2);
    }
};