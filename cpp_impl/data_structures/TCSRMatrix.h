#pragma once
#include <vector>
#include <iostream>

#include "DataStructureInterface.hpp" // If conforming

class TCSRMatrix : public DataStructureInterface
{
public:
  int num_matrix_rows;
  int num_matrix_cols;
  std::vector<int> row_offsets;  // Size = num_matrix_rows + 1. Start index in encoded_cols for each row.
  std::vector<int> encoded_cols; // Stores encoded column_index (value is +/-1).

  TCSRMatrix() : num_matrix_rows(0), num_matrix_cols(0) {}
  TCSRMatrix(const int *W_raw, int K_rows, int N_cols)
  {
    init(W_raw, K_rows, N_cols);
  }

  static int tcsr_encode_col(int col_idx, int val)
  {
    // Assumes col_idx >= 0
    if (val == 1)
      return col_idx;

    // if (val == -1)
    return ~col_idx;
  }

  static std::pair<int, int> tcsr_decode_col(int encoded_c)
  {
    if (encoded_c < 0)
    {
      return {-encoded_c - 1, -1};
    }
    return {encoded_c, 1}; // Original column, value 1
  }

  void init(const int *matrix, int rows, int cols) override
  {
    num_matrix_rows = rows;
    num_matrix_cols = cols;

    row_offsets.assign(num_matrix_rows + 1, 0);
    encoded_cols.clear();

    int current_nnz_count = 0;
    for (int r = 0; r < num_matrix_rows; ++r)
    {
      row_offsets[r] = current_nnz_count;
      for (int c = 0; c < num_matrix_cols; ++c)
      {
        int val = matrix[r * num_matrix_cols + c];
        if (val == 1 || val == -1)
        {
          encoded_cols.push_back(tcsr_encode_col(c, val));
          current_nnz_count++;
        }
      }
    }
    row_offsets[num_matrix_rows] = current_nnz_count;
  }

  std::vector<int> getVectorRepresentation(size_t expected_rows, size_t expected_cols) override
  {
    std::vector<int> M(num_matrix_rows * num_matrix_cols, 0);
    for (int r = 0; r < num_matrix_rows; ++r)
    {
      int start_offset = row_offsets[r];
      int end_offset = row_offsets[r + 1];
      for (int i = start_offset; i < end_offset; ++i)
      {
        std::pair<int, int> decoded = tcsr_decode_col(encoded_cols[i]);
        int c = decoded.first;
        int val = decoded.second;
        M[r * num_matrix_cols + c] = val;
      }
    }
    return M;
  }

  int getNumRows() const { return num_matrix_rows; }
  int getNumCols() const { return num_matrix_cols; }

  // For debugging
  void printVars() override
  {
    std::cout << "\nTCSRMatrix (" << num_matrix_rows << "x" << num_matrix_cols << "):" << std::endl;
    std::cout << "row_offsets (" << row_offsets.size() << "): ";
    for (int off : row_offsets)
      std::cout << off << " ";
    std::cout << std::endl;
    std::cout << "encoded_cols (" << encoded_cols.size() << "): ";
    for (int ec : encoded_cols)
    {
      auto dec = tcsr_decode_col(ec);
      std::cout << ec << "->(" << dec.first << "," << dec.second << ") ";
    }
    std::cout << std::endl;
  }
};