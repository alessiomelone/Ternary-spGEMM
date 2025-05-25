#pragma once
#include <vector>
#include <iostream>

// ASSUMING K DIVIDES B
template <int B>
class BlockedTCSC_interleaved
{
public:
    std::vector<int> col_start;
    std::vector<int> row_index;
    std::vector<int> excess_start;
    std::vector<int> excess_type;

    BlockedTCSC_interleaved(int *matrix, int K, int N)
    {
        // Assuming K DIVIDES B
        int num_blocks_rows = K / B;
        int column_start_val = 0;     // Renamed for clarity from original 'column_start' local variable
        int excess_start_pos_val = 0; // Renamed for clarity from original 'excess_start_pos' local variable

        for (int b = 0; b < num_blocks_rows; b++)
        {
            for (int j = 0; j < N; j++)
            {
                col_start.push_back(column_start_val);
                excess_start.push_back(excess_start_pos_val);

                std::vector<int> pos_indices;
                std::vector<int> neg_indices;

                for (int i = 0; i < B; i++)
                {
                    if (matrix[(b * B + i) * N + j] == 1)
                    {
                        pos_indices.push_back(b * B + i);
                    }
                    else if (matrix[(b * B + i) * N + j] == -1)
                    {
                        neg_indices.push_back(b * B + i);
                    }
                }

                int min_count = std::min(pos_indices.size(), neg_indices.size());
                for (int k_idx = 0; k_idx < min_count; k_idx++)
                {
                    row_index.push_back(pos_indices[k_idx]);
                    row_index.push_back(neg_indices[k_idx]);
                    column_start_val += 2; // This local variable accumulates total interleaved elements processed so far
                }

                if (pos_indices.size() > min_count)
                {
                    excess_type.push_back(1);
                    for (size_t k_idx = min_count; k_idx < pos_indices.size(); k_idx++)
                    {
                        row_index.push_back(pos_indices[k_idx]);
                        excess_start_pos_val++; // This local variable accumulates total positive excess elements processed so far
                    }
                }
                else if (neg_indices.size() > min_count)
                {
                    excess_type.push_back(-1);
                    for (size_t k_idx = min_count; k_idx < neg_indices.size(); k_idx++)
                    {
                        row_index.push_back(neg_indices[k_idx]);
                        excess_start_pos_val++; // This local variable accumulates total negative excess elements processed so far
                    }
                }
                else
                {
                    excess_type.push_back(0);
                }
            }
        }
        col_start.push_back(column_start_val);
        excess_start.push_back(excess_start_pos_val);
    }

    int getDataStructureSize() const
    {
        return sizeof(int) * (col_start.size() + row_index.size() + excess_start.size() + excess_type.size());
    }
};
