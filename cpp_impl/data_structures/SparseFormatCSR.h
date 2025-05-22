#pragma once
#include <vector>
#include <iostream>
#include "DataStructureInterface.hpp"

class SparseFormatCSR : public DataStructureInterface
{
public:
    int num_rows;
    int num_cols;
    std::vector<int> row_start_pos;
    std::vector<int> row_start_neg;
    std::vector<int> col_index_pos;
    std::vector<int> col_index_neg;

    SparseFormatCSR() : num_rows(0), num_cols(0) {}

    SparseFormatCSR(const int *matrix, int K, int N) : num_rows(K), num_cols(N)
    {
        init(matrix, K, N);
    }

    void init(const int *matrix, int rows, int cols) override
    {
        num_rows = rows;
        num_cols = cols;

        int row_start_pos_val = 0;
        int row_start_neg_val = 0;

        row_start_pos.clear();
        row_start_neg.clear();
        col_index_pos.clear();
        col_index_neg.clear();

        for (int k_idx = 0; k_idx < rows; k_idx++)
        {
            row_start_pos.push_back(row_start_pos_val);
            row_start_neg.push_back(row_start_neg_val);

            for (int n_idx = 0; n_idx < cols; n_idx++)
            {
                if (matrix[k_idx * cols + n_idx] >= 1)
                {
                    row_start_pos_val++;
                    col_index_pos.push_back(n_idx);
                }
                else if (matrix[k_idx * cols + n_idx] <= -1)
                {
                    row_start_neg_val++;
                    col_index_neg.push_back(n_idx);
                }
            }
        }

        row_start_pos.push_back(row_start_pos_val);
        row_start_neg.push_back(row_start_neg_val);
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

    int getNumRows() const override { return num_rows; }
    int getNumCols() const override { return num_cols; }

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
};