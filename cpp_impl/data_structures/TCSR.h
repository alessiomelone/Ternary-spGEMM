#pragma once
#include <vector>
#include <iostream>

class TCSR
{
public:
    std::vector<int> row_start_pos;
    std::vector<int> row_start_neg;
    std::vector<int> col_index_pos;
    std::vector<int> col_index_neg;

    TCSR(const int *matrix, int rows, int cols)
    {
        int row_start_pos_count = 0;
        int row_start_neg_count = 0;

        for (int k = 0; k < rows; k++)
        {
            row_start_pos.push_back(row_start_pos_count);
            row_start_neg.push_back(row_start_neg_count);

            for (int n = 0; n < cols; n++)
            {
                int val = matrix[k * cols + n];
                if (val == 1)
                {
                    row_start_pos_count++;
                    col_index_pos.push_back(n);
                }
                else if (val == -1)
                {
                    row_start_neg_count++;
                    col_index_neg.push_back(n);
                }
            }
        }

        row_start_pos.push_back(row_start_pos_count);
        row_start_neg.push_back(row_start_neg_count);
    }

    int getDataStructureSize() const
    {
        return sizeof(int) * (row_start_pos.size() +
                              row_start_neg.size() +
                              col_index_pos.size() +
                              col_index_neg.size());
    }
};