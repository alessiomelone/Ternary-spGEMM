#pragma once
#include <vector>
#include <iostream>

class BaseTCSC
{
public:
    std::vector<int> col_start_pos;
    std::vector<int> col_start_neg;
    std::vector<int> row_index_pos;
    std::vector<int> row_index_neg;

    BaseTCSC(const int *matrix, int rows, int cols)
    {
        int column_start_pos = 0;
        int column_start_neg = 0;

        for (int n = 0; n < cols; n++)
        {
            col_start_pos.push_back(column_start_pos);
            col_start_neg.push_back(column_start_neg);

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
        }

        col_start_pos.push_back(column_start_pos);
        col_start_neg.push_back(column_start_neg);
    }

    int getDataStructureSize() const
    {
        return sizeof(int) * (col_start_pos.size() +
                              col_start_neg.size() +
                              row_index_pos.size() +
                              row_index_neg.size());
    }
};