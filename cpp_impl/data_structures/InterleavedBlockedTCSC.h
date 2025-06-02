/* Note [Harry]: Implementation heavily inspired by Bilal's. */
#pragma once
#include <vector>
#include <iostream>

// ASSUMING K DIVIDES B
template <int B>
class InterleavedBlockedTCSC
{
public:
    std::vector<int> all_indices;     // All row indices, grouped by pattern
    std::vector<int> col_segment_ptr; // Pointers to segments within all_indices

    // matrix W is KxN
    InterleavedBlockedTCSC(const int *matrix, int rows, int cols)
    {
        int num_blocks_rows = rows / B;

        col_segment_ptr.push_back(0);

        std::vector<int> temp_pos_indices;
        std::vector<int> temp_neg_indices;

        // For each of the K/B blocks of the rows of W
        for (int block = 0; block < num_blocks_rows; block++)
        {
            int start_row = block * B;
            int end_row = start_row + B;
            // For every column of W
            for (int n = 0; n < cols; ++n)
            {
                temp_pos_indices.clear();
                temp_neg_indices.clear();

                // For every row in the block
                for (int k = start_row; k < end_row; ++k)
                {
                    int val = matrix[k * cols + n];
                    if (val == 1)
                    {
                        temp_pos_indices.push_back(k);
                    }
                    else if (val == -1)
                    {
                        temp_neg_indices.push_back(k);
                    }
                }

                size_t p_idx = 0;
                size_t n_idx = 0;
                // BASE implementation uses +,- interleaving (i.e. groups of 1 sign)
                // change group size here by chaning the + 1 to + 3 for grourps of 4
                while (p_idx < temp_pos_indices.size() && n_idx < temp_neg_indices.size())
                {
                    all_indices.push_back(temp_pos_indices[p_idx++]);
                    all_indices.push_back(temp_neg_indices[n_idx++]);
                }
                col_segment_ptr.push_back(all_indices.size());

                // Remaining positives
                while (p_idx < temp_pos_indices.size())
                {
                    all_indices.push_back(temp_pos_indices[p_idx++]);
                }
                col_segment_ptr.push_back(all_indices.size());

                // Remaining negatives
                while (n_idx < temp_neg_indices.size())
                {
                    all_indices.push_back(temp_neg_indices[n_idx++]);
                }
                col_segment_ptr.push_back(all_indices.size());
            }
        }
    }

    // matrix W is KxN
    InterleavedBlockedTCSC(const int *matrix, int rows, int cols, int UNROLL_FACTOR)
    {
        int num_blocks_rows = rows / B;

        col_segment_ptr.push_back(0);

        std::vector<int> temp_pos_indices;
        std::vector<int> temp_neg_indices;

        // For each of the K/B blocks of the rows of W
        for (int block = 0; block < num_blocks_rows; block++)
        {
            int start_row = block * B;
            int end_row = start_row + B;
            // For every column of W
            for (int n = 0; n < cols; ++n)
            {
                temp_pos_indices.clear();
                temp_neg_indices.clear();

                // For every row in the block
                for (int k = start_row; k < end_row; ++k)
                {
                    int val = matrix[k * cols + n];
                    if (val == 1)
                    {
                        temp_pos_indices.push_back(k);
                    }
                    else if (val == -1)
                    {
                        temp_neg_indices.push_back(k);
                    }
                }

                size_t p_idx = 0;
                size_t n_idx = 0;

                // UNROLLED impl uses groups of (UNROLL_FACTOR/2) signs, i.e. ++++,---- if UF=8
                // change group size here by chaning the + 1 to + 3 for grourps of 4
                while (p_idx + UNROLL_FACTOR / 2 - 1 < temp_pos_indices.size() && n_idx + UNROLL_FACTOR / 2 - 1 < temp_neg_indices.size())
                {
                    for (int i = 0; i < UNROLL_FACTOR / 2; i++)
                    {
                        // copy and add extra lines for larger or fewer groups
                        all_indices.push_back(temp_pos_indices[p_idx++]);
                    }

                    for (int i = 0; i < UNROLL_FACTOR / 2; i++)
                    {
                        // copy and add extra lines for larger or fewer groups
                        all_indices.push_back(temp_neg_indices[n_idx++]);
                    }
                }
                col_segment_ptr.push_back(all_indices.size());

                // Remaining positives
                while (p_idx < temp_pos_indices.size())
                {
                    all_indices.push_back(temp_pos_indices[p_idx++]);
                }
                col_segment_ptr.push_back(all_indices.size());

                // Remaining negatives
                while (n_idx < temp_neg_indices.size())
                {
                    all_indices.push_back(temp_neg_indices[n_idx++]);
                }
                col_segment_ptr.push_back(all_indices.size());
            }
        }
    }

    int getDataStructureSize() const
    {
        size_t size_bytes = 0;
        size_bytes += sizeof(int) * all_indices.size();
        size_bytes += sizeof(int) * col_segment_ptr.size();
        return size_bytes;
    }
};