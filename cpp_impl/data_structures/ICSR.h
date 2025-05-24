#pragma once
#include <vector>
#include <iostream>
#include <utility>
#include "DataStructureInterface.hpp"

class ICSR : public DataStructureInterface
{
public:
  int num_matrix_rows;
  int num_matrix_cols;
  std::vector<int> row_offsets;  // Size = num_matrix_rows + 1
  std::vector<int> encoded_cols; // Stores encoded column indices (±1)

  ICSR() : num_matrix_rows(0), num_matrix_cols(0) {}

  ICSR(const int *W_raw, int rows, int cols)
      : num_matrix_rows(0), num_matrix_cols(0)
  {
    init(W_raw, rows, cols);
  }

  static int icsr_encode_col(int c, int val)
  {
    return val == 1 ? c : ~c;
  }

  static std::pair<int, int> icsr_decode_col(int encoded_c)
  {
    return encoded_c < 0 ? std::make_pair(-encoded_c - 1, -1) : std::make_pair(encoded_c, 1);
  }

  void init(const int *matrix, int rows, int cols) override
  {
    num_matrix_rows = rows;
    num_matrix_cols = cols;
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

  std::vector<int> getVectorRepresentation(size_t /*expected_rows*/, size_t /*expected_cols*/) override
  {
    std::vector<int> M(num_matrix_rows * num_matrix_cols, 0);
    if (num_matrix_rows == 0 || num_matrix_cols == 0)
      return M;

    for (int r = 0; r < num_matrix_rows; ++r)
    {
      int start_idx = row_offsets[r];
      int end_idx = row_offsets[r + 1];
      for (int i = start_idx; i < end_idx; ++i)
      {
        auto [c, val] = icsr_decode_col(encoded_cols[i]);
        if (c < num_matrix_cols)
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
    std::cout << "\nICSR (" << num_matrix_rows << "x" << num_matrix_cols << "):\n";
    std::cout << "row_offsets: ";
    for (int o : row_offsets)
      std::cout << o << " ";
    std::cout << "\nencoded_cols: ";
    for (int ec : encoded_cols)
      std::cout << ec << " ";
    std::cout << "\n";
  }

  int getDataStructureSize() const
  {
    return sizeof(int) * (2 + row_offsets.size() + encoded_cols.size());
  }
};
