#pragma once
#include <vector>
#include <iostream>
#include "DataStructureInterface.hpp"

using namespace std;

class InterleavedTCSC_baraq : public DataStructureInterface
{
public:
    vector<int> col_start_pos;
    vector<int> col_start_neg;
    vector<int> row_index_pos;
    vector<int> row_index_neg;

    vector<int> row_index_interleaved;
    vector<int> col_start_interleaved;

    InterleavedTCSC_baraq() {}
    InterleavedTCSC_baraq(const int *W_raw, int K, int N)
    {
        init(W_raw, K, N);
    }

    void init(const int *matrix, int rows, int cols)
    {
        int column_start_pos = 0;
        int column_start_neg = 0;
        int column_start_interleaved = 0;

        for (int n = 0; n < cols; n++)
        {
            col_start_pos.push_back(column_start_pos);
            col_start_neg.push_back(column_start_neg);
            col_start_interleaved.push_back(column_start_interleaved);

            for (int k = 0; k < rows; k++)
            {
                int val = matrix[k * cols + n];
                if (val == 1)
                {
                    column_start_pos++;
                    row_index_pos.push_back(k);
                }
                else if (val == -1)
                {
                    column_start_neg++;
                    row_index_neg.push_back(k);
                }
            }

            // calculate number of elements to interleave
            size_t first_row_index_in_col_pos = col_start_pos[col_start_pos.size() - 1];
            size_t first_row_index_in_col_neg = col_start_neg[col_start_neg.size() - 1];
            size_t pos_elements_in_col = column_start_pos - first_row_index_in_col_pos;
            size_t neg_elements_in_col = column_start_neg - first_row_index_in_col_neg;
            size_t elements_to_interleave = min(pos_elements_in_col, neg_elements_in_col);

            // add elements to interleaved
            size_t pos_ptr = first_row_index_in_col_pos;
            size_t neg_ptr = first_row_index_in_col_neg;
            size_t i;
            for (i = 0; i < elements_to_interleave; ++i)
            {
                row_index_interleaved.push_back(row_index_pos[pos_ptr++]); // pos then neg
                row_index_interleaved.push_back(row_index_neg[neg_ptr++]);
            }

            // remove the last i elements we added to interleaved
            row_index_pos.erase(row_index_pos.begin() + first_row_index_in_col_pos, row_index_pos.begin() + first_row_index_in_col_pos + i); // from the current col_start
            row_index_neg.erase(row_index_neg.begin() + first_row_index_in_col_neg, row_index_neg.begin() + first_row_index_in_col_neg + i);
            column_start_pos -= i;
            column_start_neg -= i;
            column_start_interleaved += 2 * i;
        }

        // last element in col_start
        col_start_pos.push_back(column_start_pos);
        col_start_neg.push_back(column_start_neg);
        col_start_interleaved.push_back(column_start_interleaved);
    }

    int getDataStructureSize() const
    {
        return sizeof(int) * (col_start_pos.size() +
                              col_start_neg.size() +
                              row_index_pos.size() +
                              row_index_neg.size());
    }

    std::vector<int> getVectorRepresentation(size_t rows, size_t cols)
    {
        std::vector<int> dense(rows * cols, 0);

        const size_t stored_cols = col_start_pos.size() - 1;
        if (stored_cols != cols)
            throw std::invalid_argument("cols does not match stored matrix");

        for (size_t n = 0; n < cols; ++n)
        {
            size_t inter_start = col_start_interleaved[n];
            size_t inter_end = col_start_interleaved[n + 1];
            for (size_t idx = inter_start; idx < inter_end; idx += 2)
            {
                dense[row_index_interleaved[idx] * cols + n] = 1;
                dense[row_index_interleaved[idx + 1] * cols + n] = -1;
            }

            size_t pos_start = col_start_pos[n];
            size_t pos_end = col_start_pos[n + 1];
            for (size_t idx = pos_start; idx < pos_end; ++idx)
                dense[row_index_pos[idx] * cols + n] = 1;

            size_t neg_start = col_start_neg[n];
            size_t neg_end = col_start_neg[n + 1];
            for (size_t idx = neg_start; idx < neg_end; ++idx)
                dense[row_index_neg[idx] * cols + n] = -1;
        }

        return dense;
    }
};