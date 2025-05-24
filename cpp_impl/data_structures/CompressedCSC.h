#pragma once
#include <vector>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <random>
#include <iomanip>

#include "DataStructureInterface.hpp"
#include "compress.h"

using namespace std;

class CompressedCSC : public DataStructureInterface
{
public:
	// FOR VALUES OF 3 OR MORE
	// Each value is a byte representing 5 ternary values.
	vector<uint8_t> vals_5;
	// Contains indices in vals_5 where the next matrix column starts, inclusive.
	// Each byte is only in one column.
	vector<int> col_start_5;
	vector<int> row_index_5;

	vector<int> col_start_pos;
	vector<int> col_start_neg;
	vector<int> row_index_pos;
	vector<int> row_index_neg;

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
			col_start_5.push_back(vals_5.size());
			col_start_neg.push_back(row_index_neg.size());
			col_start_pos.push_back(row_index_pos.size());

			for (int row = 0; row < rows; ) {
				// check if the next 5 values have 3+ nonzero elements
				if (row < rows - 5) {
					int nonzeros = 0;
					for (int i = 0; i < 5; ++i)
						if (matrix[(row + i) * cols + col] != 0)
							++nonzeros;
					if (nonzeros >= 3) {
						// encode next 5 values as a block of 5
						const int *matrix_ptr = matrix + (row * cols) + col;
						uint8_t bitstring = encode5(matrix_ptr[0], matrix_ptr[1 * cols], matrix_ptr[2 * cols], matrix_ptr[3 * cols], matrix_ptr[4 * cols]);
						vals_5.push_back(bitstring);
						row_index_5.push_back(row);
						rows += 5;
						continue;
					}
				}

				// TODO: otherwise, check if the next 2 values are both nonzero

				// otherwise, encode the current value if it's nonzero
				if (matrix[row * cols + col] == 1) {
					row_index_pos.push_back(row);
				} else if (matrix[row * cols + col] == -1) {
					row_index_neg.push_back(row);
				}
				++row;
			}

			// for (int row = 0; row <= rows - 5; row += 5)
			// {
			// 	const int *matrix_ptr = matrix + (row * cols) + col;
			// 	if (matrix_ptr[0] == 0 && matrix_ptr[1 * cols] == 0 && matrix_ptr[2 * cols] == 0 && matrix_ptr[3 * cols] == 0 && matrix_ptr[4 * cols] == 0)
			// 	{
			// 		continue; // Skip bytes that are fully zero
			// 	}

			// 	uint8_t bitstring = encode5(matrix_ptr[0], matrix_ptr[1 * cols], matrix_ptr[2 * cols], matrix_ptr[3 * cols], matrix_ptr[4 * cols]);

			// 	vals_5.push_back(bitstring);
			// 	row_index_5.push_back(row);
			// }

			// cleanup for last few bits
			// int pad_values[5] = {0, 0, 0, 0, 0};
			// for (int pad_pos = 0; pad_pos < rows % 5; ++pad_pos)
			// {
			// 	int row = (rows - (rows % 5)) + pad_pos;
			// 	const int *matrix_ptr = matrix + (row * cols) + col;
			// 	pad_values[pad_pos] = *matrix_ptr;
			// }
			// uint8_t bitstring;
			// if (pad_values[0] != 0 || pad_values[1] != 0 || pad_values[2] != 0 || pad_values[3] != 0 || pad_values[4] != 0)
			// {
			// 	bitstring = encode5(pad_values[0], pad_values[1], pad_values[2], pad_values[3], pad_values[4]);
			// 	vals_5.push_back(bitstring);
			// 	row_index_5.push_back(rows - (rows % 5));
			// }
		}
		col_start_5.push_back(vals_5.size()); // end col_start_5 val
		col_start_neg.push_back(row_index_neg.size());
		col_start_pos.push_back(row_index_pos.size());
	}

	

	void printVars()
	{
		std::cout << "\nvals_5: ";
		for (uint8_t val : vals_5)
		{
			std::cout << int(val) << " ";
		}
		std::cout << std::endl;
		std::cout << "vals_5 decoded: ";
		for (auto b : vals_5)
		{
			std::cout << "[";
			for (int8_t val : decode5[b])
			{
				std::cout << int(val) << " ";
			}
			std::cout << "]";
		}
		std::cout << "\ncol_start_5: ";
		for (auto idx : col_start_5)
		{
			std::cout << idx << " ";
		}
		std::cout << "\nrow_start: ";
		for (auto idx : row_index_5)
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
			size_t start = col_start_5[c];
			// …and end (next col’s start, or end of vals_5)
			size_t end = (c + 1 < col_start_5.size()
							  ? col_start_5[c + 1]
							  : vals_5.size());

			for (size_t idx = start; idx < end; ++idx)
			{
				uint8_t code = vals_5[idx];
				int baseRow = row_index_5[idx];
				auto &dec = decode5[code]; // decoded 5-tuple

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

            // reconstruct singular +1 entries stored in row_index_pos
            size_t startPos = col_start_pos[c];
            size_t endPos = (c + 1 < col_start_pos.size() ? col_start_pos[c + 1] : row_index_pos.size());
            for (size_t idxPos = startPos; idxPos < endPos; ++idxPos)
            {
                size_t r = row_index_pos[idxPos];
                if (r < rows)
                {
                    M[r * cols + c] = 1;
                }
            }

            // reconstruct singular -1 entries stored in row_index_neg
            size_t startNeg = col_start_neg[c];
            size_t endNeg = (c + 1 < col_start_neg.size() ? col_start_neg[c + 1] : row_index_neg.size());
            for (size_t idxNeg = startNeg; idxNeg < endNeg; ++idxNeg)
            {
                size_t r = row_index_neg[idxNeg];
                if (r < rows)
                {
                    M[r * cols + c] = -1;
                }
            }
			
		}

		return M;
	}

    int getDataStructureSize() const {
        return sizeof(uint8_t) * vals_5.size() +
               sizeof(int) * (col_start_5.size() + col_start_pos.size() + col_start_neg.size()) +
               sizeof(int) * (row_index_5.size() + row_index_pos.size() + row_index_neg.size()) +
               256 * 5 * sizeof(uint8_t);
    }
};