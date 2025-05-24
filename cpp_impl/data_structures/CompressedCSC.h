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
	// Each value is a byte representing 5 ternary values.
	vector<uint8_t> vals;

	// Contains indices in vals where the next matrix column starts, inclusive.
	// Each byte is only in one column.
	vector<int> col_start;
	vector<int> row_index;

	CompressedCSC() {}
	CompressedCSC(const int *W_raw, int K, int N)
	{
		init(W_raw, K, N);
	}

	void init(const int *matrix, int rows_int, int cols_int)
	{
		int rows = static_cast<int>(rows_int);
		int cols = static_cast<int>(cols_int);
		for (int col = 0; col < cols; ++col)
		{
			col_start.push_back(vals.size());
			for (int row = 0; row <= rows - 5; row += 5)
			{
				const int *matrix_ptr = matrix + (row * cols) + col;
				if (matrix_ptr[0] == 0 && matrix_ptr[1 * cols] == 0 && matrix_ptr[2 * cols] == 0 && matrix_ptr[3 * cols] == 0 && matrix_ptr[4 * cols] == 0)
				{
					continue; // Skip bytes that are fully zero
				}

				uint8_t bitstring = encode(matrix_ptr[0], matrix_ptr[1 * cols], matrix_ptr[2 * cols], matrix_ptr[3 * cols], matrix_ptr[4 * cols]);

				vals.push_back(bitstring);
				row_index.push_back(row);
			}

			// cleanup for last few bits
			int pad_values[5] = {0, 0, 0, 0, 0};
			for (int pad_pos = 0; pad_pos < rows % 5; ++pad_pos)
			{
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
				int val = matrix[r * cols + c];
				if (val != 0)
				{
					row_index.push_back(r);
					vals.push_back(static_cast<uint8_t>(val));
					current_nnz_count++;
				}
=======
				short row = (rows - (rows % 5)) + pad_pos;
=======
				int row = (rows - (rows % 5)) + pad_pos;
>>>>>>> 67ba190 (short to int)
=======
				int row = (rows - (rows % 5)) + pad_pos;
>>>>>>> a757576 (short to int)
				const int *matrix_ptr = matrix + (row * cols) + col;
				pad_values[pad_pos] = *matrix_ptr;
			}
			uint8_t bitstring;
			if (pad_values[0] != 0 || pad_values[1] != 0 || pad_values[2] != 0 || pad_values[3] != 0 || pad_values[4] != 0)
			{
				bitstring = encode(pad_values[0], pad_values[1], pad_values[2], pad_values[3], pad_values[4]);
				vals.push_back(bitstring);
				row_index.push_back(rows - (rows % 5));
>>>>>>> 82b3402 (use shorts in CCSC & run unroll test)
			}
		}
<<<<<<< HEAD
		col_start[cols] = current_nnz_count;
=======
		col_start.push_back(vals.size()); // end col_start val
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations)
	}

	static constexpr uint8_t encode(int8_t v0, int8_t v1, int8_t v2,
									int8_t v3, int8_t v4) noexcept
	{
		return uint8_t((v0 + 1) * 81 +
					   (v1 + 1) * 27 +
					   (v2 + 1) * 9 +
					   (v3 + 1) * 3 +
					   (v4 + 1) * 1);
	}

	void printVars()
	{
		std::cout << "\nvals: ";
		for (uint8_t val : vals)
		{
			std::cout << int(val) << " ";
		}
		std::cout << std::endl;
		std::cout << "vals decoded: ";
		for (auto b : vals)
		{
			std::cout << "[";
			for (int8_t val : decodeCCSC[b])
			{
				std::cout << int(val) << " ";
			}
			std::cout << "]";
		}
		std::cout << "\ncol_start: ";
		for (auto idx : col_start)
		{
			std::cout << idx << " ";
		}
		std::cout << "\nrow_start: ";
		for (auto idx : row_index)
		{
			std::cout << idx << " ";
		}
		std::cout << std::endl;
	}

	vector<int> getVectorRepresentation(size_t rows, size_t cols)
	{
		// allocate and zero-init output
		vector<int> M(rows * cols, 0);

		// for each column
		for (size_t c = 0; c < cols; ++c)
		{
			// where this column’s bytes start…
			size_t start = col_start[c];
			// …and end (next col’s start, or end of vals)
			size_t end = (c + 1 < col_start.size()
							  ? col_start[c + 1]
							  : vals.size());

			for (size_t idx = start; idx < end; ++idx)
			{
				uint8_t code = vals[idx];
				int baseRow = row_index[idx];
				auto &dec = decodeCCSC[code]; // decoded 5-tuple

				// assign each of the up-to-5 entries, skip past-end rows
				for (int j = 0; j < 5; ++j)
				{
					size_t r = baseRow + j;
					if (r < rows)
					{
						M[r * cols + c] = dec[j];
					}
				}
			}
		}

		return M;
	}

    int getDataStructureSize() const {
        return sizeof(uint8_t) * vals.size() +
			   sizeof(int) * col_start.size() +
			   sizeof(int) * row_index.size() +
			   256 * 5 * sizeof(uint8_t);
	}
};