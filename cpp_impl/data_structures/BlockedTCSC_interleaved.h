#pragma once
#include <vector>
#include <iostream>

// ASSUMING K DIVIDES B
template <int B>
class BlockedTCSC
{
public:
    std::vector<int> col_start_pos;
    std::vector<int> col_start_neg;
    std::vector<int> row_index_pos;
    std::vector<int> row_index_neg;

    BlockedTCSC(int *matrix, int K, int N)
    {
        int num_blocks_rows = K / B;
        int column_start_pos = 0;
        int column_start_neg = 0;
        for (int b = 0; b < num_blocks_rows; b++)
        {
            for (int j = 0; j < N; j++)
            {
                col_start_pos.push_back(column_start_pos);
                col_start_neg.push_back(column_start_neg);
                for (int i = 0; i < B; i++)
                {
                    if (matrix[(b * B + i) * N + j] == 1)
                    {
                        row_index_pos.push_back(b * B + i);
                        column_start_pos++;
                    }
                    else if (matrix[(b * B + i) * N + j] == -1)
                    {
                        row_index_neg.push_back(b * B + i);
                        column_start_neg++;
                    }
                }
            }
        }
        col_start_pos.push_back(column_start_pos);
        col_start_neg.push_back(column_start_neg);
    }

    int getDataStructureSize() const
    {
        return sizeof(int) * (col_start_pos.size() + col_start_neg.size() + row_index_pos.size() + row_index_neg.size());
    }
};
