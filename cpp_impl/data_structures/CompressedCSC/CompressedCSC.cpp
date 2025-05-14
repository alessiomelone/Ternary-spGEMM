#include "CompressedCSC.h"
#include "../../SparseGEMM.h"
#include <iomanip>

template <typename T>
void printMatrix(const std::vector<T> &mat,
                 std::size_t rows,
                 std::size_t cols,
                 int fieldWidth = 2)
{
    for (std::size_t i = 0; i < rows; ++i)
    {
        for (std::size_t j = 0; j < cols; ++j)
        {
            const T &val = mat[i * cols + j];
            if (fieldWidth > 0)
                std::cout << std::setw(fieldWidth) << val;
            else
                std::cout << val;
            if (j + 1 < cols)
                std::cout << ' ';
        }
        std::cout << '\n';
    }
}

int main()
{
    int K = 6;
    int N = 6;
    int nonZero = 2;
    vector<int> W_truth = generateSparseMatrix<int>(K, N, nonZero, false);
    CompressedCSC compressed = CompressedCSC(&W_truth[0], K, N);
    printMatrix<int>(W_truth, K, N);
    printMatrix<int>(compressed, K, N);
}
