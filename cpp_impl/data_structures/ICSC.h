#pragma once
#include <vector>
#include <iostream>
#include <utility>

class ICSC
{
public:
  std::vector<int> col_offsets;
  std::vector<int> encoded_rows;

  ICSC(const int *matrix, int rows, int cols)
  {
    int current_nnz_count = 0;

    for (int c = 0; c < cols; c++)
    {
      col_offsets.push_back(current_nnz_count);

      for (int r = 0; r < rows; r++)
      {
        int val = matrix[r * cols + c];
        if (val == 1 || val == -1)
        {
          encoded_rows.push_back(icsc_encode_row(r, val));
          current_nnz_count++;
        }
      }
    }

    col_offsets.push_back(current_nnz_count);
  }

  static int icsc_encode_row(int r, int val)
  {
    return val == 1 ? r : ~r;
  }

  static std::pair<int, int> icsc_decode_row(int encoded_r)
  {
    return encoded_r < 0 ? std::make_pair(-encoded_r - 1, -1) : std::make_pair(encoded_r, 1);
  }

  int getDataStructureSize() const
  {
    return sizeof(int) * (col_offsets.size() + encoded_rows.size());
  }
};
