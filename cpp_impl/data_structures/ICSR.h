#pragma once
#include <vector>
#include <iostream>
#include <utility>

class ICSR
{
public:
  std::vector<int> row_offsets;  // Size = num_matrix_rows + 1
  std::vector<int> encoded_cols; // Stores encoded column indices (±1)

  ICSR(const int *matrix, int rows, int cols)
  {
    row_offsets.assign(rows + 1, 0);
    encoded_cols.clear();

    int current_nnz_count = 0;
    for (int r = 0; r < rows; ++r)
    {
      row_offsets[r] = current_nnz_count;
      for (int c = 0; c < cols; ++c)
      {
        int val = matrix[r * cols + c];
        if (val == 1 || val == -1)
        {
          encoded_cols.push_back(icsr_encode_col(c, val));
          current_nnz_count++;
        }
      }
    }
    row_offsets[rows] = current_nnz_count;
  }

  static int icsr_encode_col(int c, int val)
  {
    return val == 1 ? c : ~c;
  }

  static std::pair<int, int> icsr_decode_col(int encoded_c)
  {
    return encoded_c < 0 ? std::make_pair(-encoded_c - 1, -1) : std::make_pair(encoded_c, 1);
  }

  int getDataStructureSize() const
  {
    return sizeof(int) * (2 + row_offsets.size() + encoded_cols.size());
  }
};
