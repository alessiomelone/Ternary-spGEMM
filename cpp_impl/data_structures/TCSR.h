#pragma once
#include <vector>
#include <iostream>
#include <utility>
#include "DataStructureInterface.hpp"

class CompressedTCSR : public DataStructureInterface
{
public:
  int num_matrix_rows;
  int num_matrix_cols;
  std::vector<int> row_offsets;
  std::vector<int> encoded_cols;
  std::vector<int> row_pos_counts;

  CompressedTCSR() : num_matrix_rows(0), num_matrix_cols(0) {}
  CompressedTCSR(const int *W, int K, int N) { init(W, K, N); }

  static int tcsr_encode_col(int c, int v) { return v == 1 ? c : ~c; }

  static std::pair<int, int> tcsr_decode_col(int ec)
  {
    return ec < 0 ? std::make_pair(-ec - 1, -1) : std::make_pair(ec, 1);
  }

  void init(const int *matrix, int rows, int cols) override
  {
    num_matrix_rows = rows;
    num_matrix_cols = cols;
    row_offsets.assign(rows + 1, 0);
    row_pos_counts.assign(rows, 0);
    encoded_cols.clear();

    std::vector<int> pos, neg;
    int current_total = 0;

    for (int r = 0; r < rows; ++r)
    {
      row_offsets[r] = current_total;
      pos.clear();
      neg.clear();

      for (int c = 0; c < cols; ++c)
      {
        int val = matrix[r * cols + c];
        if (val == 1)
          pos.push_back(c);
        else if (val == -1)
          neg.push_back(tcsr_encode_col(c, val));
      }

      row_pos_counts[r] = pos.size();
      encoded_cols.insert(encoded_cols.end(), pos.begin(), pos.end());
      encoded_cols.insert(encoded_cols.end(), neg.begin(), neg.end());
      current_total += pos.size() + neg.size();
    }
    row_offsets[rows] = current_total;
  }

  std::vector<int> getVectorRepresentation(size_t, size_t) override
  {
    std::vector<int> M(num_matrix_rows * num_matrix_cols, 0);

    for (int r = 0; r < num_matrix_rows; ++r)
    {
      int offset = row_offsets[r];
      int pos = row_pos_counts[r];

      for (int i = 0; i < pos; ++i)
        M[r * num_matrix_cols + encoded_cols[offset + i]] = 1;

      for (int i = pos; i < row_offsets[r + 1] - offset; ++i)
      {
        auto [c, val] = tcsr_decode_col(encoded_cols[offset + i]);
        M[r * num_matrix_cols + c] = val;
      }
    }
    return M;
  }

  int getNumRows() const { return num_matrix_rows; }
  int getNumCols() const { return num_matrix_cols; }

  void printVars() override
  {
    std::cout << "\nTCSRMatrix (" << num_matrix_rows << "x" << num_matrix_cols << "):\n";
    std::cout << "row_offsets: ";
    for (int o : row_offsets)
      std::cout << o << " ";
    std::cout << "\nrow_pos_counts: ";
    for (int c : row_pos_counts)
      std::cout << c << " ";
    std::cout << "\nencoded_cols: ";
    for (int ec : encoded_cols)
      std::cout << ec << " ";
    std::cout << "\n";
  }

  int getDataStructureSize() {
    return sizeof(int) * (2 +
      row_offsets.size() +
      encoded_cols.size() +
      row_pos_counts.size()
    );
  }
};
