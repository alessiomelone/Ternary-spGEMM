#pragma once
#include <vector>
#include <iostream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <random>
#include <iomanip>

#include "DataStructureInterface.hpp"

using namespace std;

static constexpr int8_t decodeCCSC[256][5] = {
	{-1, -1, -1, -1, -1},
	{-1, -1, -1, -1, 0},
	{-1, -1, -1, -1, 1},
	{-1, -1, -1, 0, -1},
	{-1, -1, -1, 0, 0},
	{-1, -1, -1, 0, 1},
	{-1, -1, -1, 1, -1},
	{-1, -1, -1, 1, 0},
	{-1, -1, -1, 1, 1},
	{-1, -1, 0, -1, -1},
	{-1, -1, 0, -1, 0},
	{-1, -1, 0, -1, 1},
	{-1, -1, 0, 0, -1},
	{-1, -1, 0, 0, 0},
	{-1, -1, 0, 0, 1},
	{-1, -1, 0, 1, -1},
	{-1, -1, 0, 1, 0},
	{-1, -1, 0, 1, 1},
	{-1, -1, 1, -1, -1},
	{-1, -1, 1, -1, 0},
	{-1, -1, 1, -1, 1},
	{-1, -1, 1, 0, -1},
	{-1, -1, 1, 0, 0},
	{-1, -1, 1, 0, 1},
	{-1, -1, 1, 1, -1},
	{-1, -1, 1, 1, 0},
	{-1, -1, 1, 1, 1},
	{-1, 0, -1, -1, -1},
	{-1, 0, -1, -1, 0},
	{-1, 0, -1, -1, 1},
	{-1, 0, -1, 0, -1},
	{-1, 0, -1, 0, 0},
	{-1, 0, -1, 0, 1},
	{-1, 0, -1, 1, -1},
	{-1, 0, -1, 1, 0},
	{-1, 0, -1, 1, 1},
	{-1, 0, 0, -1, -1},
	{-1, 0, 0, -1, 0},
	{-1, 0, 0, -1, 1},
	{-1, 0, 0, 0, -1},
	{-1, 0, 0, 0, 0},
	{-1, 0, 0, 0, 1},
	{-1, 0, 0, 1, -1},
	{-1, 0, 0, 1, 0},
	{-1, 0, 0, 1, 1},
	{-1, 0, 1, -1, -1},
	{-1, 0, 1, -1, 0},
	{-1, 0, 1, -1, 1},
	{-1, 0, 1, 0, -1},
	{-1, 0, 1, 0, 0},
	{-1, 0, 1, 0, 1},
	{-1, 0, 1, 1, -1},
	{-1, 0, 1, 1, 0},
	{-1, 0, 1, 1, 1},
	{-1, 1, -1, -1, -1},
	{-1, 1, -1, -1, 0},
	{-1, 1, -1, -1, 1},
	{-1, 1, -1, 0, -1},
	{-1, 1, -1, 0, 0},
	{-1, 1, -1, 0, 1},
	{-1, 1, -1, 1, -1},
	{-1, 1, -1, 1, 0},
	{-1, 1, -1, 1, 1},
	{-1, 1, 0, -1, -1},
	{-1, 1, 0, -1, 0},
	{-1, 1, 0, -1, 1},
	{-1, 1, 0, 0, -1},
	{-1, 1, 0, 0, 0},
	{-1, 1, 0, 0, 1},
	{-1, 1, 0, 1, -1},
	{-1, 1, 0, 1, 0},
	{-1, 1, 0, 1, 1},
	{-1, 1, 1, -1, -1},
	{-1, 1, 1, -1, 0},
	{-1, 1, 1, -1, 1},
	{-1, 1, 1, 0, -1},
	{-1, 1, 1, 0, 0},
	{-1, 1, 1, 0, 1},
	{-1, 1, 1, 1, -1},
	{-1, 1, 1, 1, 0},
	{-1, 1, 1, 1, 1},
	{0, -1, -1, -1, -1},
	{0, -1, -1, -1, 0},
	{0, -1, -1, -1, 1},
	{0, -1, -1, 0, -1},
	{0, -1, -1, 0, 0},
	{0, -1, -1, 0, 1},
	{0, -1, -1, 1, -1},
	{0, -1, -1, 1, 0},
	{0, -1, -1, 1, 1},
	{0, -1, 0, -1, -1},
	{0, -1, 0, -1, 0},
	{0, -1, 0, -1, 1},
	{0, -1, 0, 0, -1},
	{0, -1, 0, 0, 0},
	{0, -1, 0, 0, 1},
	{0, -1, 0, 1, -1},
	{0, -1, 0, 1, 0},
	{0, -1, 0, 1, 1},
	{0, -1, 1, -1, -1},
	{0, -1, 1, -1, 0},
	{0, -1, 1, -1, 1},
	{0, -1, 1, 0, -1},
	{0, -1, 1, 0, 0},
	{0, -1, 1, 0, 1},
	{0, -1, 1, 1, -1},
	{0, -1, 1, 1, 0},
	{0, -1, 1, 1, 1},
	{0, 0, -1, -1, -1},
	{0, 0, -1, -1, 0},
	{0, 0, -1, -1, 1},
	{0, 0, -1, 0, -1},
	{0, 0, -1, 0, 0},
	{0, 0, -1, 0, 1},
	{0, 0, -1, 1, -1},
	{0, 0, -1, 1, 0},
	{0, 0, -1, 1, 1},
	{0, 0, 0, -1, -1},
	{0, 0, 0, -1, 0},
	{0, 0, 0, -1, 1},
	{0, 0, 0, 0, -1},
	{0, 0, 0, 0, 0},
	{0, 0, 0, 0, 1},
	{0, 0, 0, 1, -1},
	{0, 0, 0, 1, 0},
	{0, 0, 0, 1, 1},
	{0, 0, 1, -1, -1},
	{0, 0, 1, -1, 0},
	{0, 0, 1, -1, 1},
	{0, 0, 1, 0, -1},
	{0, 0, 1, 0, 0},
	{0, 0, 1, 0, 1},
	{0, 0, 1, 1, -1},
	{0, 0, 1, 1, 0},
	{0, 0, 1, 1, 1},
	{0, 1, -1, -1, -1},
	{0, 1, -1, -1, 0},
	{0, 1, -1, -1, 1},
	{0, 1, -1, 0, -1},
	{0, 1, -1, 0, 0},
	{0, 1, -1, 0, 1},
	{0, 1, -1, 1, -1},
	{0, 1, -1, 1, 0},
	{0, 1, -1, 1, 1},
	{0, 1, 0, -1, -1},
	{0, 1, 0, -1, 0},
	{0, 1, 0, -1, 1},
	{0, 1, 0, 0, -1},
	{0, 1, 0, 0, 0},
	{0, 1, 0, 0, 1},
	{0, 1, 0, 1, -1},
	{0, 1, 0, 1, 0},
	{0, 1, 0, 1, 1},
	{0, 1, 1, -1, -1},
	{0, 1, 1, -1, 0},
	{0, 1, 1, -1, 1},
	{0, 1, 1, 0, -1},
	{0, 1, 1, 0, 0},
	{0, 1, 1, 0, 1},
	{0, 1, 1, 1, -1},
	{0, 1, 1, 1, 0},
	{0, 1, 1, 1, 1},
	{1, -1, -1, -1, -1},
	{1, -1, -1, -1, 0},
	{1, -1, -1, -1, 1},
	{1, -1, -1, 0, -1},
	{1, -1, -1, 0, 0},
	{1, -1, -1, 0, 1},
	{1, -1, -1, 1, -1},
	{1, -1, -1, 1, 0},
	{1, -1, -1, 1, 1},
	{1, -1, 0, -1, -1},
	{1, -1, 0, -1, 0},
	{1, -1, 0, -1, 1},
	{1, -1, 0, 0, -1},
	{1, -1, 0, 0, 0},
	{1, -1, 0, 0, 1},
	{1, -1, 0, 1, -1},
	{1, -1, 0, 1, 0},
	{1, -1, 0, 1, 1},
	{1, -1, 1, -1, -1},
	{1, -1, 1, -1, 0},
	{1, -1, 1, -1, 1},
	{1, -1, 1, 0, -1},
	{1, -1, 1, 0, 0},
	{1, -1, 1, 0, 1},
	{1, -1, 1, 1, -1},
	{1, -1, 1, 1, 0},
	{1, -1, 1, 1, 1},
	{1, 0, -1, -1, -1},
	{1, 0, -1, -1, 0},
	{1, 0, -1, -1, 1},
	{1, 0, -1, 0, -1},
	{1, 0, -1, 0, 0},
	{1, 0, -1, 0, 1},
	{1, 0, -1, 1, -1},
	{1, 0, -1, 1, 0},
	{1, 0, -1, 1, 1},
	{1, 0, 0, -1, -1},
	{1, 0, 0, -1, 0},
	{1, 0, 0, -1, 1},
	{1, 0, 0, 0, -1},
	{1, 0, 0, 0, 0},
	{1, 0, 0, 0, 1},
	{1, 0, 0, 1, -1},
	{1, 0, 0, 1, 0},
	{1, 0, 0, 1, 1},
	{1, 0, 1, -1, -1},
	{1, 0, 1, -1, 0},
	{1, 0, 1, -1, 1},
	{1, 0, 1, 0, -1},
	{1, 0, 1, 0, 0},
	{1, 0, 1, 0, 1},
	{1, 0, 1, 1, -1},
	{1, 0, 1, 1, 0},
	{1, 0, 1, 1, 1},
	{1, 1, -1, -1, -1},
	{1, 1, -1, -1, 0},
	{1, 1, -1, -1, 1},
	{1, 1, -1, 0, -1},
	{1, 1, -1, 0, 0},
	{1, 1, -1, 0, 1},
	{1, 1, -1, 1, -1},
	{1, 1, -1, 1, 0},
	{1, 1, -1, 1, 1},
	{1, 1, 0, -1, -1},
	{1, 1, 0, -1, 0},
	{1, 1, 0, -1, 1},
	{1, 1, 0, 0, -1},
	{1, 1, 0, 0, 0},
	{1, 1, 0, 0, 1},
	{1, 1, 0, 1, -1},
	{1, 1, 0, 1, 0},
	{1, 1, 0, 1, 1},
	{1, 1, 1, -1, -1},
	{1, 1, 1, -1, 0},
	{1, 1, 1, -1, 1},
	{1, 1, 1, 0, -1},
	{1, 1, 1, 0, 0},
	{1, 1, 1, 0, 1},
	{1, 1, 1, 1, -1},
	{1, 1, 1, 1, 0},
	{1, 1, 1, 1, 1}};

constexpr uint8_t encode(int8_t v0, int8_t v1, int8_t v2, int8_t v3, int8_t v4) noexcept;

class CompressedCSC : public DataStructureInterface
{
public:
	int num_rows;
	int num_cols;
	std::vector<int> col_start; // Size = num_cols + 1. Start index in row_index for each column.
	std::vector<int> row_index; // Stores row indices for non-zero elements.
	std::vector<uint8_t> vals;	// Stores encoded values for non-zero elements.

	CompressedCSC() : num_rows(0), num_cols(0) {}

	CompressedCSC(const int *matrix, int rows, int cols) : num_rows(rows), num_cols(cols)
	{
		init(matrix, rows, cols);
	}

	void init(const int *matrix, int rows, int cols) override
	{
		num_rows = rows;
		num_cols = cols;

		col_start.assign(cols + 1, 0);
		row_index.clear();
		vals.clear();

		int current_nnz_count = 0;
		for (int c = 0; c < cols; ++c)
		{
			col_start[c] = current_nnz_count;
			for (int r = 0; r < rows; ++r)
			{
				int val = matrix[r * cols + c];
				if (val != 0)
				{
					row_index.push_back(r);
					vals.push_back(static_cast<uint8_t>(val));
					current_nnz_count++;
				}
			}
		}
<<<<<<< HEAD
		col_start[cols] = current_nnz_count;
=======
		col_start.push_back(vals.size()); // end col_start val
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations)
	}

	static constexpr uint8_t encode(int8_t v0, int8_t v1, int8_t v2, int8_t v3, int8_t v4) noexcept
	{
		return uint8_t((v0 + 1) * 81 +
					   (v1 + 1) * 27 +
					   (v2 + 1) * 9 +
					   (v3 + 1) * 3 +
					   (v4 + 1) * 1);
	}

	void printVars() override
	{
		std::cout << "\nCompressedCSC (" << num_rows << "x" << num_cols << "):" << std::endl;
		std::cout << "col_start: ";
		for (int val : col_start)
			std::cout << val << " ";
		std::cout << "\nrow_index: ";
		for (int val : row_index)
			std::cout << val << " ";
		std::cout << "\nvals: ";
		for (uint8_t val : vals)
			std::cout << static_cast<int>(val) << " ";
		std::cout << std::endl;
	}

	std::vector<int> getVectorRepresentation(size_t expected_rows, size_t expected_cols) override
	{
		std::vector<int> result(expected_rows * expected_cols, 0);

		for (int c = 0; c < num_cols && c < expected_cols; ++c)
		{
			int start_offset = col_start[c];
			int end_offset = col_start[c + 1];
			for (int i = start_offset; i < end_offset; ++i)
			{
				int r = row_index[i];
				if (r < expected_rows)
				{
					result[r * expected_cols + c] = static_cast<int>(vals[i]);
				}
			}
		}

		return result;
	}

	int getNumRows() const override { return num_rows; }
	int getNumCols() const override { return num_cols; }
};