#pragma once
#include <vector>
#include <iostream>
#include <numeric>
#include <algorithm>

class InterleavedTCSC
{
public:
  std::vector<int> all_indices;     // All row indices, grouped by pattern
  std::vector<int> col_segment_ptr; // Pointers to segments within all_indices

  int num_original_cols; // Number of columns

  InterleavedTCSC(const int *matrix, int rows, int cols) : num_original_cols(cols)
  {
    col_segment_ptr.push_back(0);

    std::vector<int> temp_pos_indices;
    std::vector<int> temp_neg_indices;

    for (int n = 0; n < cols; ++n)
    {
      temp_pos_indices.clear();
      temp_neg_indices.clear();

      for (int k = 0; k < rows; ++k)
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

      // change group size here by chaning the + 1 to + 3 for grourps of 4
      while (p_idx + 3 < temp_pos_indices.size() && n_idx + 3 < temp_neg_indices.size())
      {
        // copy and add extra lines for larger or fewer groups
        all_indices.push_back(temp_pos_indices[p_idx++]);
        all_indices.push_back(temp_pos_indices[p_idx++]);
        all_indices.push_back(temp_pos_indices[p_idx++]);
        all_indices.push_back(temp_pos_indices[p_idx++]);

        // copy and add extra lines for larger or fewer groups
        all_indices.push_back(temp_neg_indices[n_idx++]);
        all_indices.push_back(temp_neg_indices[n_idx++]);
        all_indices.push_back(temp_neg_indices[n_idx++]);
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

  int getDataStructureSize() const
  {
    size_t size_bytes = 0;
    size_bytes += sizeof(int) * all_indices.size();
    size_bytes += sizeof(int) * col_segment_ptr.size();
    size_bytes += sizeof(int);
    return size_bytes;
  }
};