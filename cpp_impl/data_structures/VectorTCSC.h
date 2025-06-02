#pragma once
#include <vector>
#include <iostream>

class VectorTCSC
{
public:
    std::vector<int> row_index_pos;
    std::vector<int> row_index_neg;
    std::vector<int> cap_every_four;

    VectorTCSC(const int *matrix, int rows, int cols)
    {
        int cap = 0;
        // Calculate max number of either +1 or -1 in the next 4 columns
        for (int n = 0; n < cols; n++)
        {
            if (n % 4 == 0)
            {
                int max = 0;
                for (int col = n; col < cols; col++)
                {
                    int pos = 0, min = 0;
                    for (int row = 0; row < rows; row++) 
                    {
                        int val = matrix[row * cols + col];
                        if (val == 1)
                            pos++;
                        else if (val == -1)
                            min++;
                    }
                    int local_max = pos >= min ? pos : min;
                    if (local_max > max)
                        max = local_max;
                }
                cap = max;
                while (cap % 4) cap++; // make it div by 4
                cap_every_four.push_back(cap);
            }

            int pos_added = 0, minus_added = 0;
            for (int k = 0; k < rows; k++)
            {
                int val = matrix[k * cols + n];
                if (val == 1)
                {
                    pos_added++;
                    row_index_pos.push_back(k);
                }
                else if (val == -1)
                {
                    minus_added++;
                    row_index_neg.push_back(k);
                }
            }
            while (cap > pos_added++)
            {
                row_index_pos.push_back(-1);
            }
            while (minus_added++ < cap)
            {
                row_index_neg.push_back(-1);
            }
        }
    }

    int getDataStructureSize() const
    {
        return sizeof(int) * (cap_every_four.size() +
                              row_index_pos.size() +
                              row_index_neg.size());
    }
};