#pragma once
#include <vector>
#include <iostream>
#include <utility> // For std::pair
#include <numeric> // For std::iota, std::fill

#include "DataStructureInterface.hpp" // Assuming you want it to conform

class TCSCMatrix : public DataStructureInterface
{
public:
  int num_matrix_rows;
  int num_matrix_cols;
  std::vector<int> col_offsets;  // Size: num_matrix_cols + 1. Start index in encoded_rows for each col.
  std::vector<int> encoded_rows; // Stores encoded row_index (value is +/-1).

  TCSCMatrix() : num_matrix_rows(0), num_matrix_cols(0) {}
  TCSCMatrix(const int *W_raw, int K_rows, int N_cols)
  {
    init(W_raw, K_rows, N_cols);
  }

  // Encodes a plain row index and its sign (+1 or -1)
  static int tcsc_encode_row(int row_idx, int val)
  {
    // Assumes row_idx >= 0
    if (val == 1)
      return row_idx;

    return ~row_idx;
  }

  // Decodes an encoded row index back to original row and its sign
  static std::pair<int, int> tcsc_decode_row(int encoded_r)
  {
    if (encoded_r < 0)
    {
      return {-encoded_r - 1, -1}; // Original row, value -1
    }
    return {encoded_r, 1}; // Original row, value 1
  }

  void init(const int *matrix, int rows, int cols) override
  {
    num_matrix_rows = rows;
    num_matrix_cols = cols;

    col_offsets.assign(num_matrix_cols + 1, 0);
    encoded_rows.clear();

    int current_nnz_count = 0;
    for (int c = 0; c < num_matrix_cols; ++c)
    { // Iterate column by column
      col_offsets[c] = current_nnz_count;
      for (int r = 0; r < num_matrix_rows; ++r)
      {
        int val = matrix[r * num_matrix_cols + c];
        if (val == 1 || val == -1)
        {
          encoded_rows.push_back(tcsc_encode_row(r, val));
          current_nnz_count++;
        }
      }
    }
    col_offsets[num_matrix_cols] = current_nnz_count;
  }

  std::vector<int> getVectorRepresentation(size_t expected_rows, size_t expected_cols) override
  {
    std::vector<int> M(num_matrix_rows * num_matrix_cols, 0);
    if (num_matrix_rows == 0 || num_matrix_cols == 0)
      return M;

    for (int c = 0; c < num_matrix_cols; ++c)
    { // Iterate columns
      int start_offset = col_offsets[c];
      int end_offset = col_offsets[c + 1];
      for (int i = start_offset; i < end_offset; ++i)
      {
        std::pair<int, int> decoded = tcsc_decode_row(encoded_rows[i]);
        int r = decoded.first;
        int val = decoded.second;
        if (r >= 0 && r < num_matrix_rows)
        {
          M[r * num_matrix_cols + c] = val; // M[r][c] = val
        }
      }
    }
    return M;
  }

  void printVars()
  {
    std::cout << "\nTCSCMatrix (" << num_matrix_rows << "x" << num_matrix_cols << "):" << std::endl;
    std::cout << "col_offsets (" << col_offsets.size() << "): ";
    for (int off : col_offsets)
      std::cout << off << " ";
    std::cout << std::endl;
    std::cout << "encoded_rows (" << encoded_rows.size() << "): ";
    for (int er : encoded_rows)
    {
      auto dec = tcsc_decode_row(er);
      std::cout << er << "->(" << dec.first << "," << dec.second << ") ";
    }
    std::cout << std::endl;
  }
};