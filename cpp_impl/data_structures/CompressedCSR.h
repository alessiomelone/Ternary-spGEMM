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

class CompressedCSR : public DataStructureInterface
{
public:
    // Each value is a byte representing 5 ternary values.
    vector<uint8_t> vals;

    // Contains indices in vals where the next matrix column starts, inclusive.
    // Each byte is only in one column.
    vector<int> row_start;
    vector<int> col_index;

    CompressedCSR() {}
    CompressedCSR(const int *W_raw, int K, int N)
    {
        init(W_raw, K, N);
    }

    void init(const int *matrix, int rows, int cols)
    {
        for (int row = 0; row < rows; ++row)
        {
            row_start.push_back(vals.size());
            for (int col = 0; col <= cols - 5; col += 5)
            {
                const int *matrix_ptr = matrix + (row * cols) + col;
                if (matrix_ptr[0] == 0 && matrix_ptr[1] == 0 && matrix_ptr[2] == 0 && matrix_ptr[3] == 0 && matrix_ptr[4] == 0)
                {
                    continue; // Skip bytes that are fully zero
                }

                uint8_t bitstring = encode5(matrix_ptr[0], matrix_ptr[1], matrix_ptr[2], matrix_ptr[3], matrix_ptr[4]);
                vals.push_back(bitstring);
                col_index.push_back(col);
            }

            // cleanup for last few bits
            int pad_values[5] = {0, 0, 0, 0, 0};
            for (int pad_pos = 0; pad_pos < cols % 5; ++pad_pos)
            {
                int col = (cols - (cols % 5)) + pad_pos;
                const int *matrix_ptr = matrix + (row * cols) + col;
                pad_values[pad_pos] = *matrix_ptr;
            }
            uint8_t bitstring;
            if (pad_values[0] != 0 || pad_values[1] != 0 || pad_values[2] != 0 || pad_values[3] != 0 || pad_values[4] != 0)
            {
                bitstring = encode5(pad_values[0], pad_values[1], pad_values[2], pad_values[3], pad_values[4]);
                vals.push_back(bitstring);
                col_index.push_back(cols - (cols % 5));
            }
        }
        row_start.push_back(vals.size()); // end col_start val
    }

    /*
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
            for (int8_t val : decode5[b])
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
    */
    vector<int> getVectorRepresentation(size_t rows, size_t cols)
    {
        vector<int> M(rows * cols, 0);

        for (size_t r = 0; r < rows; ++r)
        {
            size_t start = row_start[r];   // first packed byte in this row
            size_t end = row_start[r + 1]; // one-past-last byte in this row

            for (size_t idx = start; idx < end; ++idx)
            {
                const uint8_t code = vals[idx];
                const int baseCol = col_index[idx];   // column of the first value in this byte
                const int8_t *dec = decode5[code]; // 5 decoded ternary values

                for (int j = 0; j < 5; ++j)
                {
                    size_t c = baseCol + j;
                    if (c < cols) // guard for last partial byte
                        M[r * cols + c] = dec[j];
                }
            }
        }
        return M;
    }

    int getDataStructureSize() const
    {
        return sizeof(uint8_t) * vals.size() +
               sizeof(int) * row_start.size() +
               sizeof(int) * col_index.size() +
               256 * 5 * sizeof(uint8_t);
    }
};