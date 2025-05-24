#pragma once
#include <vector>
#include <iostream>
#include <utility>
#include "DataStructureInterface.hpp"

class ICSC : public DataStructureInterface
{
public:
  int num_matrix_rows;
  int num_matrix_cols;
  std::vector<int> col_offsets;
  std::vector<int> encoded_rows;

  ICSC() : num_matrix_rows(0), num_matrix_cols(0) {}

  ICSC(const int *W_raw, int rows, int cols)
      : num_matrix_rows(0), num_matrix_cols(0)
  {
    init(W_raw, rows, cols);
  }

  static int icsc_encode_row(int r, int val)
  {
    return val == 1 ? r : ~r;
  }

  static std::pair<int, int> icsc_decode_row(int encoded_r)
  {
    return encoded_r < 0 ? std::make_pair(-encoded_r - 1, -1) : std::make_pair(encoded_r, 1);
  }

  void init(const int *matrix, int rows, int cols) override
  {
    num_matrix_rows = rows;
    num_matrix_cols = cols;
    col_offsets.assign(cols + 1, 0);
    encoded_rows.clear();

    int current_nnz_count = 0;
    for (int c = 0; c < cols; ++c)
    {
      col_offsets[c] = current_nnz_count;
      for (int r = 0; r < rows; ++r)
      {
        int val = matrix[r * cols + c];
        if (val == 1 || val == -1)
        {
          encoded_rows.push_back(icsc_encode_row(r, val));
          current_nnz_count++;
        }
      }
    }
    col_offsets[cols] = current_nnz_count;
  }

  std::vector<int> getVectorRepresentation(size_t expected_rows, size_t expected_cols) override
  {
    std::vector<int> M(num_matrix_rows * num_matrix_cols, 0);
    if (num_matrix_rows == 0 || num_matrix_cols == 0)
      return M;

    for (int c = 0; c < num_matrix_cols; ++c)
    {
      int start_idx = col_offsets[c];
      int end_idx = col_offsets[c + 1];
      for (int i = start_idx; i < end_idx; ++i)
      {
        auto [r, val] = icsc_decode_row(encoded_rows[i]);
        if (r < num_matrix_rows)
        {
          M[r * num_matrix_cols + c] = val;
        }
      }
    }
    return M;
  }

  int getNumRows() const { return num_matrix_rows; }
  int getNumCols() const { return num_matrix_cols; }

  void printVars() override
  {
    std::cout << "\nICSC (" << num_matrix_rows << "x" << num_matrix_cols << "):\n";
    std::cout << "col_offsets: ";
    for (int o : col_offsets)
      std::cout << o << " ";
    std::cout << "\nencoded_rows: ";
    for (int er : encoded_rows)
      std::cout << er << " ";
    std::cout << "\n";
  }

  int getDataStructureSize() const
  {
    return sizeof(int) * (2 +
                          col_offsets.size() +
                          encoded_rows.size());
  }
};
