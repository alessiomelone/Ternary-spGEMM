#pragma once
#include <vector>
#include <iostream>
#include <algorithm>

class InterleavedTCSC
{
public:
    std::vector<int> col_starts;
    std::vector<int> interleaved_row_indices;
    std::vector<int> interleaved_signs;

    InterleavedTCSC(const int *matrix, int rows, int cols)
    {

        std::vector<int> temp_row_index_pos;
        std::vector<int> temp_row_index_neg;

        int nz_count = 0;

        for (int n = 0; n < cols; ++n)
        {
            col_starts.push_back(nz_count);

            temp_row_index_pos.clear();
            temp_row_index_neg.clear();

            for (int k = 0; k < rows; ++k)
            {
                int val = matrix[k * cols + n];
                if (val == 1)
                {
                    temp_row_index_pos.push_back(k);
                }
                else if (val == -1)
                {
                    temp_row_index_neg.push_back(k);
                }
            }

            size_t p_idx = 0;
            size_t n_idx = 0;
            size_t len_pos = temp_row_index_pos.size();
            size_t len_neg = temp_row_index_neg.size();

            while (p_idx + 1 < len_pos && n_idx + 1 < len_neg)
            {
                interleaved_row_indices.push_back(temp_row_index_pos[p_idx++]);
                interleaved_signs.push_back(1);
                interleaved_row_indices.push_back(temp_row_index_pos[p_idx++]);
                interleaved_signs.push_back(1);

                interleaved_row_indices.push_back(temp_row_index_neg[n_idx++]);
                interleaved_signs.push_back(-1);
                interleaved_row_indices.push_back(temp_row_index_neg[n_idx++]);
                interleaved_signs.push_back(-1);
                nz_count += 4;
            }

            while (p_idx < len_pos)
            {
                interleaved_row_indices.push_back(temp_row_index_pos[p_idx++]);
                interleaved_signs.push_back(1);
                nz_count++;
            }

            while (n_idx < len_neg)
            {
                interleaved_row_indices.push_back(temp_row_index_neg[n_idx++]);
                interleaved_signs.push_back(-1);
                nz_count++;
            }
        }
        col_starts.push_back(nz_count);
    }

    int getDataStructureSize() const
    {
        return sizeof(int) * (col_starts.size() +
                              interleaved_row_indices.size() +
                              interleaved_signs.size());
    }
};