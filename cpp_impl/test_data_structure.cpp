#include "data_structures/CompressedCSC.h"
#include "sparseUtils.h"
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

bool vecsEqual(const vector<int> &v1, const vector<int> &v2)
{
    if (v1.size() != v2.size())
        return false;

    for (size_t i = 0; i < v1.size(); ++i)
    {
        if (v1[i] != v2[i])
            return false;
    }

    return true;
}

template <typename T>
bool test(int K, int N, int nonZero, int seed, bool verbose = false)
{
    static_assert(std::is_base_of<DataStructureInterface, T>::value,
                  "T must inherit from DataStructureInterface");

    auto W_truth = generateSparseMatrix<int>(K, N, nonZero, false, seed);

    if (verbose)
    {
        printMatrix<int>(W_truth, K, N, 2);
        printf("------------------------------\n");
    }

    T matrix_representation;
    matrix_representation.init(&W_truth[0], K, N);
    auto M = matrix_representation.getVectorRepresentation(K, N);

    if (verbose)
    {
        printMatrix<int>(M, K, N, 2);
        printf(M == W_truth ? "pass\n" : "fail\n");
    }

    return M == W_truth;
}

template <typename T>
bool testMany(int K, int N, int nonZero, int variants, bool verbose = false)
{
    static_assert(std::is_base_of<DataStructureInterface, T>::value,
                  "T must inherit from DataStructureInterface");

    // const int K = 10, N = 10, nonZero = 2, variants = 1'000;
    bool mismatch_flag = false;
    for (int k = 1; k < K; ++k)
    {
        for (int n = 1; n < N; ++n)
        {

            if (verbose)
            {
                printf("Testing k=%d, n=%d, nonZero=%d\n", k, n, nonZero);
            }

            for (int seed = 0; seed < variants; ++seed)
            {
                if (!test<T>(k, n, nonZero, seed))
                {
                    printf("Mismatch at k=%d, n=%d, seed=%d\n", k, n, seed);
                    mismatch_flag = true;
                    return false;
                }
            }
        }
    }
    if (!mismatch_flag)
        printf("All vectors match!\n");
    return true;
}

template <typename T>
bool testRequired(int variants, bool verbose = false)
{
    static_assert(std::is_base_of<DataStructureInterface, T>::value,
                  "T must inherit from DataStructureInterface");

    int K[] = {512, 1024, 2048, 4096, 2048, 4096, 8192, 16384};
    int N[] = {2048, 4096, 8192, 16384, 512, 1024, 2048, 4096};
    int nonZeros[] = {2, 4, 8, 16};

    bool mismatch_flag = false;
    for (int i = 0; i < 8; ++i)
    {
        for (int nonZero : nonZeros)
        {
            if (verbose)
            {
                printf("Testing k=%d, n=%d, nonZero=%d\n", K[i], N[i], nonZero);
            }

            for (int seed = 0; seed < variants; ++seed)
            {

                if (!test<T>(K[i], N[i], nonZero, seed))
                {
                    printf("Mismatch at k=%d, n=%d, nonZero=%d, seed=%d\n", K[i], N[i], nonZero, seed);
                    mismatch_flag = true;
                    break;
                }
            }
        }
    }
    if (!mismatch_flag)
        printf("All vectors match!\n");
    return 0;
}

int main()
{
    // test<CompressedCSC>(1, 1, 2, 1, true);

    // Sizes you'll probably use debugging
    printf("(1/2)\n");
    testMany<CompressedCSC>(40, 40, 2, 20, false);

    // Required sizes
    printf("(2/2)\n");
    testRequired<CompressedCSC>(10, true);
}
