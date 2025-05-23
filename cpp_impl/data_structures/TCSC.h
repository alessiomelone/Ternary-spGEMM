#pragma once
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include "DataStructureInterface.hpp"

class CompressedTCSC : public DataStructureInterface
{
public:
  int num_matrix_rows;
  int num_matrix_cols;
  std::vector<int> col_offsets;
  std::vector<int> encoded_rows;
  std::vector<int> col_pos_counts;

  CompressedTCSC() : num_matrix_rows(0), num_matrix_cols(0) {}
  CompressedTCSC(const int *W_raw, int K, int N) { init(W_raw, K, N); }

  static int tcsc_encode_row(int r, int val)
  {
    return val == 1 ? r : ~r;
  }

  static std::pair<int, int> tcsc_decode_row(int encoded_r)
  {
    return encoded_r < 0 ? std::make_pair(-encoded_r - 1, -1) : std::make_pair(encoded_r, 1);
  }

  void init(const int *matrix, int rows, int cols) override
  {
    num_matrix_rows = rows;
    num_matrix_cols = cols;
    col_offsets.assign(cols + 1, 0);
    col_pos_counts.assign(cols, 0);
    encoded_rows.clear();

    std::vector<int> pos, neg;
    int current_total = 0;

    for (int c = 0; c < cols; ++c)
    {
      col_offsets[c] = current_total;
      pos.clear();
      neg.clear();

      for (int r = 0; r < rows; ++r)
      {
        int val = matrix[r * cols + c];
        if (val == 1)
          pos.push_back(r);
        else if (val == -1)
          neg.push_back(tcsc_encode_row(r, val));
      }

      col_pos_counts[c] = pos.size();
      encoded_rows.insert(encoded_rows.end(), pos.begin(), pos.end());
      encoded_rows.insert(encoded_rows.end(), neg.begin(), neg.end());
      current_total += pos.size() + neg.size();
    }
    col_offsets[cols] = current_total;
  }

  std::vector<int> getVectorRepresentation(size_t, size_t) override
  {
    std::vector<int> M(num_matrix_rows * num_matrix_cols, 0);

    for (int c = 0; c < num_matrix_cols; ++c)
    {
      int offset = col_offsets[c];
      int pos = col_pos_counts[c];

      for (int i = 0; i < pos; ++i)
        M[encoded_rows[offset + i] * num_matrix_cols + c] = 1;

      for (int i = pos; i < col_offsets[c + 1] - offset; ++i)
      {
        auto [r, val] = tcsc_decode_row(encoded_rows[offset + i]);
        M[r * num_matrix_cols + c] = val;
      }
    }
    return M;
  }

  int getNumRows() const { return num_matrix_rows; }
  int getNumCols() const { return num_matrix_cols; }

  void printVars() override
  {
    std::cout << "\nTCSCMatrix (" << num_matrix_rows << "x" << num_matrix_cols << "):\n";
    std::cout << "col_offsets: ";
    for (int o : col_offsets)
      std::cout << o << " ";
    std::cout << "\ncol_pos_counts: ";
    for (int c : col_pos_counts)
      std::cout << c << " ";
    std::cout << "\nencoded_rows: ";
    for (int er : encoded_rows)
      std::cout << er << " ";
    std::cout << "\n";
  }

  int getDataStructureSize() {
    return sizeof(int) * (2 +
           col_offsets.size() +
           encoded_rows.size() +
           col_pos_counts.size()
    );
  }
};
